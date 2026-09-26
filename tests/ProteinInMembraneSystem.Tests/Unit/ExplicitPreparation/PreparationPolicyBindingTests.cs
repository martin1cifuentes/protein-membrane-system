using System.Collections.Immutable;
using System.Security.Cryptography;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class PreparationPolicyBindingTests
{
    [Fact]
    public void Exact_preparation_scope_rejects_changed_source_chemistry_leaflet_or_conditions()
    {
        var source = new StructuralSource("source", SourceRouteKind.Rcsb, "identified", "source.pdb",
            new string('a', 64), "6QWR", "first model");
        var intended = new IntendedProteinModel("intended", source, 0, null,
            ImmutableArray.Create(new ChainSelection("A", "A")),
            ImmutableArray<PartnerSelection>.Empty, ImmutableArray<AlternateLocationChoice>.Empty);
        var prepared = new AssessedPreparedProtein("prepared", "revision", intended,
            new MolecularArtifact("molecule", "prepared.pdb", new string('c', 64), "graph.json",
                null, null, 3205, null, new string('b', 64)), "chemical", ImmutableArray.Create(
                new ResidueVariantChoice(new ResidueAddress(0, "A", 108, "", "A"), "HID", "review")),
            ImmutableArray<PreparationChangeProposal>.Empty,
            new SourceToResultCorrespondence("source", "molecule", ImmutableArray<AtomCorrespondence>.Empty, true),
            ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<ScientificFinding>.Empty,
            ImmutableArray<string>.Empty, ChemicalStatePolicyVersion: "1",
            StructuralAssessmentPolicyId: "structural", StructuralAssessmentPolicyVersion: "1");
        var conditions = FixedStudyConditions.Initial;
        var upper = new LeafletComposition(LeafletSide.Upper,
            ImmutableArray.Create(new LipidFraction("DMPC", 1)));
        var lower = new LeafletComposition(LeafletSide.Lower,
            ImmutableArray.Create(new LipidFraction("DMPC", 1)));
        var model = new MembraneModel("model", upper, lower, conditions, "controlled bilayer");
        var membrane = new AssessedMembraneModel("assessed-membrane", "revision", model,
            ImmutableArray<MolecularRepresentation>.Empty, ImmutableArray<ScientificEvidence>.Empty,
            ImmutableArray<string>.Empty, PolicyId: "membrane", PolicyVersion: "1");
        var oriented = prepared.Molecule with { Id = "oriented" };
        var proposal = new PlacementProposal("proposal", prepared.Id, model.Id,
            ProteinTopologyKind.MembraneSpanning, oriented, 0, 20, 14.6,
            PlacementPhysicalSide.Both, "extracellular upper", ImmutableArray<string>.Empty,
            ImmutableArray<ScientificEvidence>.Empty, ImmutableArray<string>.Empty);
        var revision = new StudyRevision("revision", 1, intended, model, proposal.Id, conditions);
        var scope = new PreparationPolicyScope(source.Sha256, 0, null, intended.Chains,
            ImmutableArray<string>.Empty, prepared.Molecule.TopologySha256!, "chemical", "1",
            PreparationPolicyFingerprint.ComputeResidueVariants(prepared.ResidueVariants), upper, lower,
            conditions, ProteinTopologyKind.MembraneSpanning,
            "structural", "1", "membrane", "1");
        Assert.True(ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation.PolicyScopeMatches(scope, revision, prepared, membrane, proposal));
        Assert.False(ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation.PolicyScopeMatches(null, revision, prepared, membrane, proposal));
        var invalidScopes = new[]
        {
            scope with { SourceCoordinateSha256 = new string('0', 64) },
            scope with { SourceModelIndex = 1 },
            scope with { BiologicalAssemblyId = "other" },
            scope with { ChainCopies = ImmutableArray.Create(new ChainSelection("B", "A")) },
            scope with { ChainCopies = scope.ChainCopies.Add(scope.ChainCopies[0]) },
            scope with { RetainedPartnerSourceIds = ImmutableArray.Create("partner") },
            scope with { PreparedBondGraphSha256 = new string('0', 64) },
            scope with { ChemicalStatePolicyId = "other" },
            scope with { ChemicalStatePolicyVersion = "2" },
            scope with { ProteinStructuralPolicyId = "other" },
            scope with { ProteinStructuralPolicyVersion = "2" },
            scope with { MembraneSupportPolicyId = "other" },
            scope with { MembraneSupportPolicyVersion = "2" },
            scope with { ResidueVariantsSha256 = new string('0', 64) },
            scope with { Upper = upper with { Fractions = ImmutableArray.Create(new LipidFraction("POPC", 1)) } },
            scope with { Lower = lower with { PhysicalSide = LeafletSide.Upper } },
            scope with { Conditions = conditions with { OptionalTemperatureKelvin = 310 } },
            scope with { TopologyKind = ProteinTopologyKind.OneSurfaceAssociated }
        };
        foreach (var changed in invalidScopes)
            Assert.False(ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation.PolicyScopeMatches(changed, revision, prepared,
                membrane, proposal));
        var duplicateVariants = prepared with { ResidueVariants = prepared.ResidueVariants.Add(
            prepared.ResidueVariants[0] with { Variant = "HIE" }) };
        Assert.False(ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation.PolicyScopeMatches(scope, revision, duplicateVariants,
            membrane, proposal));
    }

    [Fact]
    public void Exact_6qwr_source_and_declared_chemical_choices_have_a_reproducible_variant_digest()
    {
        var directory = new DirectoryInfo(AppContext.BaseDirectory);
        while (directory is not null && !File.Exists(Path.Combine(directory.FullName,
                   "config", "policies", "source-assets", "6QWR.pdb")))
            directory = directory.Parent;
        var source = Path.Combine(directory?.FullName ?? string.Empty,
            "config", "policies", "source-assets", "6QWR.pdb");
        Assert.True(File.Exists(source), "The exact identified source asset is required for this catalogue digest proof.");
        Assert.Equal("f4c1503a60321c0cfe513e8211e43e20e2199aa5f71f22429e0e0ae8c97b0779",
            Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(source))).ToLowerInvariant());
        var variants = Read6QwrVariantChoices(source);
        Assert.Equal(35, variants.Length);
        Assert.Single(variants.Where(item => item.Variant == "HID" && item.Residue.Residue == 108));
        Assert.Equal("6b124e55f47c23fc57d7059d88da3adbcc80453fa334176059b8f5e037932951",
            PreparationPolicyFingerprint.ComputeResidueVariants(variants));
        Assert.Equal(PreparationPolicyFingerprint.ComputeResidueVariants(variants),
            PreparationPolicyFingerprint.ComputeResidueVariants(variants.Reverse().Select(item =>
                item with { DecisionId = "different review identity" }).ToImmutableArray()));
        Assert.Equal(string.Empty, PreparationPolicyFingerprint.ComputeResidueVariants(
            variants.Add(variants[0] with { Variant = "ASH" })));
    }

    internal static ImmutableArray<ResidueVariantChoice> Read6QwrVariantChoices(string source)
    {
        var residues = File.ReadLines(source).Where(line => line.StartsWith("ATOM", StringComparison.Ordinal))
            .Select(line => new
            {
                Name = line.Substring(17, 3).Trim(),
                Chain = line.Substring(21, 1),
                Number = int.Parse(line.Substring(22, 4), System.Globalization.CultureInfo.InvariantCulture),
                InsertionCode = line.Substring(26, 1).Trim()
            })
            .DistinctBy(item => (item.Chain, item.Number, item.InsertionCode)).ToArray();
        Assert.Equal(211, residues.Length);
        var defaults = new Dictionary<string, string>(StringComparer.Ordinal)
        { ["ASP"] = "ASP", ["GLU"] = "GLU", ["CYS"] = "CYS", ["LYS"] = "LYS" };
        return residues.Select(residue =>
        {
            var variant = residue.Name == "HIS" && residue.Number == 108 ? "HID" :
                defaults.GetValueOrDefault(residue.Name);
            return (residue, variant);
        }).Where(item => item.variant is not null)
            .Select(item => new ResidueVariantChoice(new ResidueAddress(0, item.residue.Chain,
                item.residue.Number, item.residue.InsertionCode, "A"), item.variant!, null))
            .ToImmutableArray();
    }
}
