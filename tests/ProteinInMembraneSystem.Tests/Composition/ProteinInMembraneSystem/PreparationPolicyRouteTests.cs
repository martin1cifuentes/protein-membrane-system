using System.Collections.Immutable;
using System.Reflection;
using System.Text.Json;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using ConstructionOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class PreparationPolicyRouteTests
{
    [Fact]
    public void Native_popc_mechanism_requires_its_own_exact_installed_patch_and_policy_species()
    {
        using var fixture = new ConstructionFixture();
        var dmpc = fixture.Policy.Construction;
        var popcPath = Path.Combine(fixture.Directory, "POPC.pdb");
        File.WriteAllText(popcPath, "controlled installed POPC patch bytes");
        var popc = dmpc with
        {
            LipidTypeArgument = "POPC", NativePatchPath = popcPath,
            NativePatchSha256 = ConstructionFixture.Hash(popcPath)
        };
        var installed = new ConstructionProviderInstallation(dmpc.ProviderVersion,
            fixture.NativePatchPath, fixture.NativePatchSha,
            ImmutableArray.Create(new NativePatchInstallation("POPC", popcPath, popc.NativePatchSha256)));
        var providerMatches = typeof(ProductRoot).GetMethod("ConstructionProviderMatches",
            BindingFlags.Static | BindingFlags.NonPublic)!;
        var validPolicy = typeof(ConstructionOwner).GetMethod("ValidConstructionPolicy",
            BindingFlags.Static | BindingFlags.NonPublic)!;
        bool Matches(ConstructionPolicy policy, ConstructionProviderInstallation observed) =>
            (bool)providerMatches.Invoke(null, [policy, observed])!;
        bool Valid(ConstructionPolicy policy) => (bool)validPolicy.Invoke(null, [policy])!;

        Assert.True(Valid(dmpc));
        Assert.DoesNotContain("nativePatchMode", JsonSerializer.Serialize(dmpc,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)), StringComparison.Ordinal);
        Assert.True(Valid(popc));
        Assert.True(Matches(dmpc, installed));
        Assert.True(Matches(popc, installed));
        Assert.False(Matches(popc, installed with { AdditionalPatches = ImmutableArray<NativePatchInstallation>.Empty }));
        Assert.False(Matches(popc, installed with { AdditionalPatches = ImmutableArray.Create(
            new NativePatchInstallation("POPC", popcPath, new string('0', 64))) }));
        Assert.False(Matches(popc, installed with { AdditionalPatches = installed.AdditionalPatches.Add(
            new NativePatchInstallation("POPC", popcPath, popc.NativePatchSha256)) }));
        Assert.False(Valid(popc with { LipidTypeArgument = "POPE" }));
        Assert.False(Matches(popc with { ProviderVersion = "other" }, installed));

        var derivedPath = Path.Combine(fixture.Directory, "POPC-63x63.pdb");
        File.WriteAllText(derivedPath, "controlled derived POPC patch bytes");
        var custom = popc with
        {
            NativePatchMode = "popc-62-109-deletion",
            NativePatchPath = derivedPath,
            NativePatchSha256 = ConstructionFixture.Hash(derivedPath),
            NativeSourcePatchPath = popcPath,
            NativeSourcePatchSha256 = popc.NativePatchSha256,
            RemovedNativeLipidResidueIds = ImmutableArray.Create("62", "109")
        };
        Assert.True(Valid(custom));
        Assert.Contains("popc-62-109-deletion", JsonSerializer.Serialize(custom,
            new JsonSerializerOptions(JsonSerializerDefaults.Web)), StringComparison.Ordinal);
        Assert.True(Matches(custom, installed));
        Assert.False(Valid(custom with { RemovedNativeLipidResidueIds = ImmutableArray.Create("109") }));
        Assert.False(Valid(custom with { NativePatchMode = "other" }));
        Assert.False(Matches(custom with { NativeSourcePatchSha256 = new string('0', 64) }, installed));
        Assert.False(Matches(custom with { NativeSourcePatchPath = derivedPath }, installed));
    }

    [Fact]
    public void Installed_slice4_catalogue_offers_only_its_exact_6qwr_dmpc_preparation_scope()
    {
        var root = RepositoryRoot();
        var cataloguePath = Path.Combine(root, "config/policies/protein-membrane-slice4.json");
        var sourcePath = Path.Combine(root, "config/policies/source-assets/6QWR.pdb");
        using var fixture = new ConstructionFixture();
        using var http = new HttpClient();
        using var catalogueDocument = JsonDocument.Parse(File.ReadAllText(cataloguePath));
        var construction = catalogueDocument.RootElement.GetProperty("preparationPolicies")[0]
            .GetProperty("construction");
        var patchPath = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(cataloguePath)!,
            construction.GetProperty("nativePatchPath").GetString()!));
        var installed = new ConstructionProviderInstallation(
            construction.GetProperty("providerVersion").GetString()!, patchPath,
            construction.GetProperty("nativePatchSha256").GetString()!);
        var product = new ProductRoot(new ScientificWorkerExchange("/missing/python", root),
            new ExternalSourceExchange(http), Path.Combine(fixture.Directory, "root-workspace"),
            cataloguePath, "/missing/ppm", () => installed);
        Assert.Null(Private<string?>(product, "_catalogueIssue"));
        var catalogue = Private<object?>(product, "_catalogue");
        Assert.NotNull(catalogue);
        var policies = (ImmutableArray<ApplicablePreparationPolicy>)catalogue.GetType()
            .GetProperty("PreparationPolicies")!.GetValue(catalogue)!;
        var policy = Assert.Single(policies);
        var lipids = (ImmutableArray<MolecularRepresentation>)catalogue.GetType()
            .GetProperty("Lipids")!.GetValue(catalogue)!;
        var dmpc = Assert.Single(lipids.Where(item => item.SpeciesId == "DMPC"));
        var source = fixture.Protein.Intended.Source with
        { CoordinatePath = sourcePath, Sha256 = ConstructionFixture.Hash(sourcePath) };
        var intended = fixture.Protein.Intended with { Source = source };
        var variants = PreparationPolicyBindingTests.Read6QwrVariantChoices(sourcePath);
        var protein = fixture.Protein with
        {
            Intended = intended,
            Molecule = fixture.Protein.Molecule with
            { TopologySha256 = policy.Scope!.PreparedBondGraphSha256 },
            ChemicalStatePolicyId = policy.Scope.ChemicalStatePolicyId,
            ChemicalStatePolicyVersion = policy.Scope.ChemicalStatePolicyVersion,
            StructuralAssessmentPolicyId = policy.Scope.ProteinStructuralPolicyId,
            StructuralAssessmentPolicyVersion = policy.Scope.ProteinStructuralPolicyVersion,
            ResidueVariants = variants
        };
        var membrane = fixture.Membrane with
        { SpeciesRepresentations = ImmutableArray.Create(dmpc),
            PolicyId = policy.Scope!.MembraneSupportPolicyId,
            PolicyVersion = policy.Scope.MembraneSupportPolicyVersion };
        var proposal = fixture.Placement.Proposal with
        { PreparedProteinId = protein.Id, MembraneModelId = membrane.Intended.Id };
        var placement = fixture.Placement with { Proposal = proposal };
        var revision = fixture.Revision with
        { IntendedProtein = intended, Membrane = membrane.Intended,
            AdoptedPlacementProposalId = proposal.Id };

        Assert.True(ConstructionOwner.PolicyScopeMatches(policy.Scope, revision, protein,
            membrane, proposal), "The installed exact source, chemical choices, graph and leaflets must match.");
        var selector = typeof(ProductRoot).GetMethod("SelectPreparationPolicyLocked",
            BindingFlags.Instance | BindingFlags.NonPublic)!;
        ApplicablePreparationPolicy? Selected(StudyRevision selectedRevision,
            AssessedPreparedProtein selectedProtein, AssessedMembraneModel selectedMembrane,
            PlacementProposal selectedProposal) =>
            (ApplicablePreparationPolicy?)selector.Invoke(product,
                [selectedRevision, selectedProtein, selectedProposal, selectedMembrane]);
        Assert.Equal(policy.Id, Selected(revision, protein, membrane, proposal)?.Id);

        SetPrivate(product, "_study", revision);
        SetPrivate(product, "_protein", protein);
        SetPrivate(product, "_membrane", membrane);
        SetPrivate(product, "_placement", placement);
        Assert.True(StartAction(product).Enabled, StartAction(product).Reason);

        var wrongProtein = protein with { Intended = intended with
            { Source = source with { Sha256 = new string('0', 64) } } };
        Assert.Null(Selected(revision, wrongProtein, membrane, proposal));
        SetPrivate(product, "_protein", wrongProtein);
        Assert.False(StartAction(product).Enabled);
        SetPrivate(product, "_protein", protein);
        Assert.Null(Selected(revision, protein with { Intended = intended with { ModelIndex = 1 } },
            membrane, proposal));
        Assert.Null(Selected(revision, protein with { Intended = intended with
            { BiologicalAssemblyId = "other" } }, membrane, proposal));
        Assert.Null(Selected(revision, protein with { Intended = intended with
            { Chains = ImmutableArray.Create(new ChainSelection("A", "B")) } }, membrane, proposal));
        Assert.Null(Selected(revision, protein with { Intended = intended with
            { Partners = ImmutableArray.Create(new PartnerSelection("partner", true, "retained")) } },
            membrane, proposal));
        Assert.Null(Selected(revision, protein with { Molecule = protein.Molecule with
            { TopologySha256 = new string('0', 64) } }, membrane, proposal));
        Assert.Null(Selected(revision, protein with { ResidueVariants = variants.SetItem(0,
            variants[0] with { Variant = "other" }) }, membrane, proposal));
        Assert.Null(Selected(revision, protein with { ChemicalStatePolicyVersion = "other" },
            membrane, proposal));
        Assert.Null(Selected(revision, protein with { StructuralAssessmentPolicyId = "other" },
            membrane, proposal));
        Assert.Null(Selected(revision, protein with { StructuralAssessmentPolicyVersion = "other" },
            membrane, proposal));
        Assert.Null(Selected(revision, protein,
            membrane with { PolicyId = "other" }, proposal));
        Assert.Null(Selected(revision, protein,
            membrane with { PolicyVersion = "other" }, proposal));
        Assert.Null(Selected(revision with { Conditions = revision.Conditions with
            { OptionalTemperatureKelvin = 310 } }, protein, membrane, proposal));
        Assert.Null(Selected(revision, protein,
            membrane with { Intended = membrane.Intended with { Lower =
                membrane.Intended.Lower with { Fractions = ImmutableArray.Create(
                    new LipidFraction("POPC", 1)) } } }, proposal));
        Assert.Null(Selected(revision, protein, membrane,
            proposal with { PreparedProteinId = "other" }));
        Assert.Null(Selected(revision, protein, membrane,
            proposal with { TopologyKind = ProteinTopologyKind.OneSurfaceAssociated }));

        var preparationPolicies = catalogue.GetType().GetProperty("PreparationPolicies")!;
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { MaximumMinimizationIterations = 0 }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { FinalUnrestrainedRmsForceTargetKjMolNm = 11 }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { Construction = policy.Construction with { MaximumConstructionSeconds = 0 } }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { Construction = policy.Construction with { MinimumPaddingNanometers = 0 } }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { Construction = policy.Construction with
                { NativePatchSha256 = new string('0', 64) } }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { ForceFieldFiles = policy.ForceFieldFiles.SetItem(0,
                policy.ForceFieldFiles[0] with { Sha256 = new string('0', 64) }) }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, ImmutableArray.Create(policy with
            { Water = policy.Water with { CoordinateTemplateSha256 = new string('0', 64) } }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        preparationPolicies.SetValue(catalogue, policies.Add(policy with { Id = "equally-applicable" }));
        Assert.Null(Selected(revision, protein, membrane, proposal));
        Assert.False(StartAction(product).Enabled);

    }

    private static AvailableAction StartAction(ProductRoot product) => product.Snapshot().Actions
        .Single(item => item.Kind == ActorActionKind.StartPreparation);

    private static T Private<T>(ProductRoot product, string name) =>
        (T)typeof(ProductRoot).GetField(name, BindingFlags.Instance | BindingFlags.NonPublic)!
            .GetValue(product)!;

    private static void SetPrivate(ProductRoot product, string name, object value) =>
        typeof(ProductRoot).GetField(name, BindingFlags.Instance | BindingFlags.NonPublic)!
            .SetValue(product, value);

    private static string RepositoryRoot()
    {
        var directory = new DirectoryInfo(AppContext.BaseDirectory);
        while (directory is not null && !File.Exists(Path.Combine(directory.FullName,
                   "config/policies/protein-membrane-slice4.json")))
            directory = directory.Parent;
        return Assert.IsType<DirectoryInfo>(directory).FullName;
    }
}
