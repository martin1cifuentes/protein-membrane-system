using System.Collections.Immutable;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using MembraneOwner = ProteinInMembrane.Host.ProteinInMembraneSystem.MembraneModelAssessment.MembraneModelAssessment;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class MembraneModelAssessmentOwnerTests
{
    [Fact]
    public async Task Zero_fraction_species_is_absent_from_policy_and_symmetric_leaflet_comparison()
    {
        var worker = new MembraneWorkerStub(request => Observed(request,
            Observation("POPC", 134)));
        var owner = new MembraneOwner(worker);
        var chosen = Model(
            ImmutableArray.Create(new LipidFraction("POPC", 1), new LipidFraction("UNSELECTED", 0)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));

        var result = await owner.AssessAsync(Revision(chosen), chosen, Catalogue(), Policy(),
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.True(result.Established, result.Reason);
        Assert.Equal("model-1", result.Value!.Intended.Id);
        Assert.Equal("revision-1", result.Value.StudyRevisionId);
        Assert.Equal(["POPC"], result.Value.SpeciesRepresentations.Select(item => item.SpeciesId));
        Assert.Single(worker.Requests);
        Assert.Equal(["POPC"], worker.Requests[0].Payload.SpeciesRepresentations.Select(item => item.SpeciesId));
        Assert.Contains("policy-popc", result.Value.Evidence.Single().Applicability);
        Assert.Contains("Membrane-local", result.Value.Evidence.Single().Uncertainty);
        Assert.Equal(("policy-popc", "1"), (result.Value.PolicyId, result.Value.PolicyVersion));
        Assert.Equal(["POPC"], worker.Requests[0].Payload.Upper.Fractions.Select(item => item.SpeciesId));
    }

    [Fact]
    public async Task Mixed_asymmetric_leaflets_remain_assessable_when_a_zero_fraction_species_is_absent()
    {
        var worker = new MembraneWorkerStub(request => Observed(request,
            new MembraneAssessmentObservations(
                ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "POPC-chemistry", 134,
                    134, true, ImmutableArray<string>.Empty),
                    new SpeciesTemplateObservation("DOPC", "DOPC-chemistry", 138,
                        138, true, ImmutableArray<string>.Empty)),
                ImmutableArray<string>.Empty, true)));
        var owner = new MembraneOwner(worker);
        var chosen = Model(
            ImmutableArray.Create(new LipidFraction("POPC", 0.6), new LipidFraction("DOPC", 0.4),
                new LipidFraction("UNSELECTED", 0)),
            ImmutableArray.Create(new LipidFraction("DOPC", 1)));
        var policy = Policy() with
        {
            CoveredSpeciesIds = ImmutableArray.Create("POPC", "DOPC"),
            AllowsMixedLeaflets = true,
            AllowsAsymmetricLeaflets = true
        };

        var result = await owner.AssessAsync(Revision(chosen), chosen, Catalogue(), policy,
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.True(result.Established, result.Reason);
        Assert.Equal(["POPC", "DOPC"], result.Value!.SpeciesRepresentations.Select(item => item.SpeciesId));
        Assert.Single(worker.Requests);
        Assert.Equal(["POPC", "DOPC"], worker.Requests[0].Payload.Upper.Fractions.Select(item => item.SpeciesId));
        Assert.Equal(["DOPC"], worker.Requests[0].Payload.Lower.Fractions.Select(item => item.SpeciesId));
    }

    [Fact]
    public async Task Same_model_id_with_changed_purpose_or_leaflet_is_not_the_adopted_revision_choice()
    {
        var worker = new MembraneWorkerStub(request => Observed(request, Observation("POPC", 134)));
        var owner = new MembraneOwner(worker);
        var adopted = Model(ImmutableArray.Create(new LipidFraction("POPC", 1)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var changedPurpose = adopted with { ScientificPurpose = "A different scientific study premise." };
        var changedLeaflet = adopted with
        {
            Lower = new LeafletComposition(LeafletSide.Lower,
                ImmutableArray.Create(new LipidFraction("DOPC", 1)))
        };
        var duplicateAdopted = adopted with
        {
            Upper = new LeafletComposition(LeafletSide.Upper,
                ImmutableArray.Create(new LipidFraction("POPC", 0.5),
                    new LipidFraction("POPC", 0.5)))
        };
        var coherentReplacement = adopted with
        {
            Upper = new LeafletComposition(LeafletSide.Upper,
                ImmutableArray.Create(new LipidFraction("POPC", 0.5),
                    new LipidFraction("DOPC", 0.5)))
        };

        foreach (var replacement in new[] { changedPurpose, changedLeaflet })
        {
            var result = await owner.AssessAsync(Revision(adopted), replacement,
                Catalogue(), Policy(), "/tmp/membrane-unit", TestContext.Current.CancellationToken);
            Assert.False(result.Established);
            Assert.Contains("study revision", result.Reason);
        }
        var duplicateCounterexample = await owner.AssessAsync(Revision(duplicateAdopted),
            coherentReplacement, Catalogue(), Policy(), "/tmp/membrane-unit",
            TestContext.Current.CancellationToken);
        Assert.False(duplicateCounterexample.Established);
        Assert.Contains("study revision", duplicateCounterexample.Reason);
        Assert.Empty(worker.Requests);
    }

    [Fact]
    public async Task Incoherent_or_uncovered_choice_does_not_invoke_scientific_worker()
    {
        var worker = new MembraneWorkerStub(request => Observed(request, Observation("POPC", 134)));
        var owner = new MembraneOwner(worker);
        var duplicate = Model(ImmutableArray.Create(new LipidFraction("POPC", 0.5),
            new LipidFraction("POPC", 0.5)), ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var unsupported = Model(ImmutableArray.Create(new LipidFraction("OTHER", 1)),
            ImmutableArray.Create(new LipidFraction("OTHER", 1)));
        var mixed = Model(ImmutableArray.Create(new LipidFraction("POPC", 0.5),
            new LipidFraction("DOPC", 0.5)), ImmutableArray.Create(new LipidFraction("POPC", 1)));

        var duplicateResult = await owner.AssessAsync(Revision(duplicate), duplicate, Catalogue(), Policy(),
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);
        var unsupportedResult = await owner.AssessAsync(Revision(unsupported), unsupported, Catalogue(), Policy(),
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);
        var mixedResult = await owner.AssessAsync(Revision(mixed), mixed, Catalogue(),
            Policy() with { CoveredSpeciesIds = ImmutableArray.Create("POPC", "DOPC") },
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.False(duplicateResult.Established);
        Assert.Contains("unambiguous", duplicateResult.Reason);
        Assert.False(unsupportedResult.Established);
        Assert.Contains("representation", unsupportedResult.Reason);
        Assert.Contains("OTHER", unsupportedResult.Reason);
        Assert.False(mixedResult.Established);
        Assert.Contains("mixture", mixedResult.Reason);
        Assert.Empty(worker.Requests);
    }

    [Theory]
    [InlineData("negative fraction", "upper", "negative")]
    [InlineData("nonfinite fraction", "upper", "finite")]
    [InlineData("nonunit total", "upper", "sum")]
    [InlineData("wrong upper physical side", "upper", "physical side")]
    [InlineData("wrong lower physical side", "lower", "physical side")]
    [InlineData("absent policy", "policy", "selected species")]
    [InlineData("disallowed asymmetry", "asymmetry", "policy")]
    public async Task Invalid_choice_or_policy_identifies_the_affected_condition_without_worker_call(
        string condition, string expectedSubject, string expectedIssue)
    {
        var worker = new MembraneWorkerStub(request => Observed(request, Observation("POPC", 134)));
        var owner = new MembraneOwner(worker);
        var chosen = Model(ImmutableArray.Create(new LipidFraction("POPC", 1)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        MembraneSupportPolicy? policy = Policy();
        switch (condition)
        {
            case "negative fraction":
                chosen = chosen with { Upper = new LeafletComposition(LeafletSide.Upper,
                    ImmutableArray.Create(new LipidFraction("POPC", -0.1), new LipidFraction("DOPC", 1.1))) };
                break;
            case "nonfinite fraction":
                chosen = chosen with { Upper = new LeafletComposition(LeafletSide.Upper,
                    ImmutableArray.Create(new LipidFraction("POPC", double.NaN))) };
                break;
            case "nonunit total":
                chosen = chosen with { Upper = new LeafletComposition(LeafletSide.Upper,
                    ImmutableArray.Create(new LipidFraction("POPC", 0.8))) };
                break;
            case "wrong upper physical side":
                chosen = chosen with { Upper = chosen.Upper with { PhysicalSide = LeafletSide.Lower } };
                break;
            case "wrong lower physical side":
                chosen = chosen with { Lower = chosen.Lower with { PhysicalSide = LeafletSide.Upper } };
                break;
            case "absent policy":
                policy = null;
                break;
            case "disallowed asymmetry":
                chosen = chosen with { Lower = new LeafletComposition(LeafletSide.Lower,
                    ImmutableArray.Create(new LipidFraction("DOPC", 1))) };
                policy = policy with { CoveredSpeciesIds = ImmutableArray.Create("POPC", "DOPC") };
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(condition));
        }

        var result = await owner.AssessAsync(Revision(chosen), chosen, Catalogue(), policy,
            "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains(expectedSubject, result.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Contains(expectedIssue, result.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Empty(worker.Requests);
    }

    [Fact]
    public async Task Each_physical_leaflet_needs_a_positive_lipid_fraction_even_when_sterol_is_parameterizable()
    {
        var worker = new MembraneWorkerStub(request => Observed(request, Observation("CHL1", 74)));
        var owner = new MembraneOwner(worker);
        var sterolOnly = ImmutableArray.Create(new LipidFraction("CHL1", 1), new LipidFraction("POPC", 0));
        var mixed = ImmutableArray.Create(new LipidFraction("POPC", 0.6), new LipidFraction("CHL1", 0.4));
        var policy = Policy() with
        {
            CoveredSpeciesIds = ImmutableArray.Create("POPC", "CHL1"),
            AllowsMixedLeaflets = true,
            AllowsAsymmetricLeaflets = true
        };

        foreach (var (side, chosen) in new[]
                 {
                     ("upper", Model(sterolOnly, mixed)),
                     ("lower", Model(mixed, sterolOnly))
                 })
        {
            var result = await owner.AssessAsync(Revision(chosen), chosen, Catalogue(), policy,
                "/tmp/membrane-unit", TestContext.Current.CancellationToken);

            Assert.False(result.Established);
            Assert.Contains(side, result.Reason, StringComparison.OrdinalIgnoreCase);
            Assert.Contains("lipid", result.Reason, StringComparison.OrdinalIgnoreCase);
        }
        Assert.Empty(worker.Requests);
    }

    [Fact]
    public async Task Uncorrelated_or_incomplete_worker_observation_cannot_establish_a_membrane()
    {
        var chosen = Model(ImmutableArray.Create(new LipidFraction("POPC", 1)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var valid = Observation("POPC", 134);
        var outcomes = new (string Label, Func<ScientificWorkRequest<MembraneAssessmentPayload>,
            WorkerResult<MembraneAssessmentObservations>> Reply)[]
        {
            ("wrong request", request => Observed(request, valid) with { RequestId = "another-request" }),
            ("wrong revision", request => Observed(request, valid) with { StudyRevisionId = "another-revision" }),
            ("absent observations", request => Observed(request, valid) with { Observations = null }),
            ("wrong chemistry", request => Observed(request, Observation("POPC", 134) with
                { Species = ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "different", 134, 134,
                    true, ImmutableArray<string>.Empty)) }) ),
            ("bond mismatch", request => Observed(request, Observation("POPC", 134) with
                { Species = ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "POPC-chemistry", 134, 134,
                    false, ImmutableArray<string>.Empty)) }) ),
            ("unobserved combination", request => Observed(request, valid with
                { CombinedParameterizationObserved = null }) ),
            ("combination warning", request => Observed(request, valid with
                { CombinationWarnings = ImmutableArray.Create("incompatible templates") }) ),
        };

        foreach (var (label, reply) in outcomes)
        {
            var worker = new MembraneWorkerStub(reply);
            var result = await new MembraneOwner(worker).AssessAsync(Revision(chosen), chosen,
                Catalogue(), Policy(), "/tmp/membrane-unit", TestContext.Current.CancellationToken);
            Assert.False(result.Established, label);
            Assert.Single(worker.Requests);
        }
    }

    [Fact]
    public async Task Exact_species_failure_and_combined_parameter_failure_identify_the_failed_aspect()
    {
        var chosen = Model(ImmutableArray.Create(new LipidFraction("POPC", 1)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var stereoWorker = new MembraneWorkerStub(request => Observed(request,
            Observation("POPC", 134) with
            {
                Species = ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "POPC-chemistry",
                    134, 134, false, ImmutableArray.Create("Tetrahedral geometry differs.")))
            }));
        var stereoResult = await new MembraneOwner(stereoWorker).AssessAsync(Revision(chosen), chosen,
            Catalogue(), Policy(), "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        var combinationWorker = new MembraneWorkerStub(request => Observed(request,
            Observation("POPC", 134) with
            {
                CombinedParameterizationObserved = false,
                CombinationWarnings = ImmutableArray.Create("Combined force-field parameterization failed.")
            }));
        var combinationResult = await new MembraneOwner(combinationWorker).AssessAsync(Revision(chosen),
            chosen, Catalogue(), Policy(), "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.False(stereoResult.Established);
        Assert.Contains("POPC", stereoResult.Reason);
        Assert.Contains("Tetrahedral geometry", stereoResult.Reason);
        Assert.False(combinationResult.Established);
        Assert.Contains("combination", combinationResult.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("parameterization failed", combinationResult.Reason);
    }

    [Theory]
    [InlineData("wrong species", "DOPC", "POPC")]
    [InlineData("unequal atom counts", "POPC", "atom counts")]
    [InlineData("worker warning", "POPC", "Unresolved stereochemistry")]
    public async Task Mismatched_worker_species_atoms_or_warnings_cannot_establish_support(
        string mismatch, string expectedSubject, string expectedIssue)
    {
        var chosen = Model(ImmutableArray.Create(new LipidFraction("POPC", 1)),
            ImmutableArray.Create(new LipidFraction("POPC", 1)));
        var worker = new MembraneWorkerStub(request => Observed(request, mismatch switch
        {
            "wrong species" => Observation("DOPC", 138),
            "unequal atom counts" => Observation("POPC", 134) with
            {
                Species = ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "POPC-chemistry",
                    134, 133, true, ImmutableArray<string>.Empty))
            },
            "worker warning" => Observation("POPC", 134) with
            {
                Species = ImmutableArray.Create(new SpeciesTemplateObservation("POPC", "POPC-chemistry",
                    134, 134, true, ImmutableArray.Create("Unresolved stereochemistry.")))
            },
            _ => throw new ArgumentOutOfRangeException(nameof(mismatch))
        }));

        var result = await new MembraneOwner(worker).AssessAsync(Revision(chosen), chosen,
            Catalogue(), Policy(), "/tmp/membrane-unit", TestContext.Current.CancellationToken);

        Assert.False(result.Established);
        Assert.Contains(expectedSubject, result.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Contains(expectedIssue, result.Reason, StringComparison.OrdinalIgnoreCase);
        Assert.Single(worker.Requests);
    }

    private static MembraneModel Model(ImmutableArray<LipidFraction> upper, ImmutableArray<LipidFraction> lower) =>
        new("model-1", new LeafletComposition(LeafletSide.Upper, upper),
            new LeafletComposition(LeafletSide.Lower, lower), FixedStudyConditions.Initial,
            "A deliberately simplified reference membrane for this study.");

    private static StudyRevision Revision(MembraneModel model) =>
        new("revision-1", 2, null, model, null, FixedStudyConditions.Initial);

    private static IReadOnlyDictionary<string, MolecularRepresentation> Catalogue() =>
        new Dictionary<string, MolecularRepresentation>(StringComparer.Ordinal)
        {
            ["POPC"] = Representation("POPC", 134),
            ["DOPC"] = Representation("DOPC", 138),
            ["CHL1"] = Representation("CHL1", 74, "sterol"),
        };

    private static MolecularRepresentation Representation(string species, int atomCount,
        string category = "lipid") => new(
        species, species + "-chemistry", category, "/identified/lipid21.xml", "a".PadLeft(64, 'a'),
        "/identified/" + species + ".pdb", "b".PadLeft(64, 'b'), atomCount, 0,
        68, 1200, ImmutableArray.Create(1, 20), "Lipid21", "8.6.0",
        ImmutableArray.Create("An intended molecular template, not an achieved bilayer."),
        ImmutableArray.Create(new MolecularStereoCheck("tetrahedral",
            ImmutableArray.Create("O21", "C1", "C3", "HS"), "negative")));

    private static MembraneSupportPolicy Policy() => new("policy-popc", "1", ImmutableArray.Create("POPC"),
        false, false, true, ImmutableArray.Create("identified local membrane-policy evidence"),
        ImmutableArray.Create("Membrane-local representation only."));

    private static MembraneAssessmentObservations Observation(string species, int atoms) => new(
        ImmutableArray.Create(new SpeciesTemplateObservation(species, species + "-chemistry", atoms,
            atoms, true, ImmutableArray<string>.Empty)), ImmutableArray<string>.Empty, true);

    private static WorkerResult<MembraneAssessmentObservations> Observed(
        ScientificWorkRequest<MembraneAssessmentPayload> request, MembraneAssessmentObservations observations) =>
        new(request.RequestId, request.Payload.StudyRevisionId, null, null, WorkerResultStanding.Observed,
            ImmutableArray<WorkerArtifact>.Empty, observations, new ProviderIdentity("controlled OpenMM facts", "8.6.0"),
            null, null);

    private sealed class MembraneWorkerStub(
        Func<ScientificWorkRequest<MembraneAssessmentPayload>, WorkerResult<MembraneAssessmentObservations>> reply)
        : IMembraneModelAssessmentWork
    {
        public List<ScientificWorkRequest<MembraneAssessmentPayload>> Requests { get; } = [];

        public Task<WorkerResult<MembraneAssessmentObservations>> AssessMembraneAsync(
            ScientificWorkRequest<MembraneAssessmentPayload> request, CancellationToken cancellationToken)
        {
            Requests.Add(request);
            return Task.FromResult(reply(request));
        }
    }
}
