using System.Net;
using System.Security.Cryptography;
using System.Text;
using ProteinInMembrane.Host;
using ProteinInMembrane.Host.ProteinInMembraneSystem;
using Xunit;

namespace ProteinInMembraneSystem.Tests;

public sealed class ExactSourceReferenceTests
{
    [Fact]
    public async Task Exact_rcsb_reference_survives_discovery_outage_and_preserves_download_identity()
    {
        var coordinates = Encoding.ASCII.GetBytes("data_exact_source\n#\n");
        var requested = new List<Uri>();
        using var http = new HttpClient(new RespondingHandler(request =>
        {
            requested.Add(request.RequestUri!);
            return request.RequestUri!.Host == "files.rcsb.org"
                ? Bytes(coordinates)
                : new HttpResponseMessage(HttpStatusCode.ServiceUnavailable);
        }));
        var exchange = new ExternalSourceExchange(http);
        var discovery = await exchange.SearchAsync("example protein", CancellationToken.None);
        Assert.Empty(discovery.Candidates);
        Assert.Equal(2, discovery.UnavailableRoutes.Length);

        var candidate = Assert.IsType<CandidateSourceRecord>(
            ExternalSourceExchange.IdentifyExactReference(SourceRouteKind.Rcsb, "1abc"));
        using var directory = new TemporaryDirectory();
        var source = await exchange.RetrieveAsync(candidate, directory.Path, CancellationToken.None);

        Assert.Equal("rcsb:1ABC", source.Id);
        Assert.Equal(SourceRouteKind.Rcsb, source.Kind);
        Assert.Equal("1ABC", source.Accession);
        Assert.Equal(coordinates, await File.ReadAllBytesAsync(source.CoordinatePath,
            TestContext.Current.CancellationToken));
        Assert.Equal(Convert.ToHexString(SHA256.HashData(coordinates)).ToLowerInvariant(), source.Sha256);
        Assert.Contains("Researcher-supplied exact RCSB PDB entry reference", source.Provenance);
        Assert.Contains(requested, uri => uri.AbsoluteUri == "https://files.rcsb.org/download/1ABC.cif");
    }

    [Theory]
    [InlineData(SourceRouteKind.Rcsb, "1ABC/other")]
    [InlineData(SourceRouteKind.Rcsb, "ABC")]
    [InlineData(SourceRouteKind.AlphaFold, "P12345")]
    [InlineData(SourceRouteKind.AlphaFold, "AF-P12345-F0")]
    public void Malformed_or_nonexact_references_are_rejected_before_network_access(
        SourceRouteKind kind, string identifier)
    {
        var candidate = ExternalSourceExchange.IdentifyExactReference(kind, identifier);
        Assert.Null(candidate);
    }

    [Fact]
    public async Task Exact_alphafold_model_selects_only_the_matching_prediction_and_preserves_evidence_scope()
    {
        var coordinates = Encoding.ASCII.GetBytes("data_prediction\n#\n");
        var metadata = """
            [
              {"modelEntityId":"AF-P12345-F2","cifUrl":"https://alphafold.ebi.ac.uk/files/second.cif"},
              {"modelEntityId":"AF-P12345-F1","sequenceStart":4,"sequenceEnd":93,
               "cifUrl":"https://alphafold.ebi.ac.uk/files/first.cif"}
            ]
            """;
        var requested = new List<Uri>();
        using var http = new HttpClient(new RespondingHandler(request =>
        {
            requested.Add(request.RequestUri!);
            return request.RequestUri!.AbsolutePath switch
            {
                "/api/prediction/P12345" => Text(metadata),
                "/files/first.cif" => Bytes(coordinates),
                _ => new HttpResponseMessage(HttpStatusCode.NotFound)
            };
        }));
        var candidate = Assert.IsType<CandidateSourceRecord>(
            ExternalSourceExchange.IdentifyExactReference(SourceRouteKind.AlphaFold, "af-p12345-f1"));
        using var directory = new TemporaryDirectory();
        var source = await new ExternalSourceExchange(http).RetrieveAsync(candidate, directory.Path, CancellationToken.None);

        Assert.Equal("alphafold:AF-P12345-F1", source.Id);
        Assert.Equal("AF-P12345-F1", source.Prediction?.RecordId);
        Assert.Equal(source.Sha256, source.Prediction?.CoordinateSha256);
        Assert.Equal(4, source.Prediction?.SequenceStart);
        Assert.Equal(93, source.Prediction?.SequenceEnd);
        Assert.Equal(PaeAcquisitionStanding.Unavailable, source.Prediction?.PaeStanding);
        Assert.Contains("did not supply a PAE", source.Prediction?.PaeReason);
        Assert.DoesNotContain(requested, uri => uri.AbsolutePath == "/files/second.cif");
    }

    [Fact]
    public async Task Duplicate_matching_alphafold_records_are_refused_before_coordinate_download()
    {
        var metadata = """
            [
              {"modelEntityId":"AF-P12345-F1","cifUrl":"https://alphafold.ebi.ac.uk/files/one.cif"},
              {"modelEntityId":"AF-P12345-F1","cifUrl":"https://alphafold.ebi.ac.uk/files/two.cif"}
            ]
            """;
        var requested = new List<Uri>();
        using var http = new HttpClient(new RespondingHandler(request =>
        {
            requested.Add(request.RequestUri!);
            return Text(metadata);
        }));
        var candidate = Assert.IsType<CandidateSourceRecord>(
            ExternalSourceExchange.IdentifyExactReference(SourceRouteKind.AlphaFold, "AF-P12345-F1"));
        using var directory = new TemporaryDirectory();
        await Assert.ThrowsAsync<InvalidDataException>(() =>
            new ExternalSourceExchange(http).RetrieveAsync(candidate, directory.Path, CancellationToken.None));
        Assert.Single(requested);
        Assert.Equal("/api/prediction/P12345", requested[0].AbsolutePath);
    }

    [Fact]
    public async Task Missing_exact_alphafold_model_is_refused_before_coordinate_download()
    {
        var metadata = """
            [{"modelEntityId":"AF-P12345-F2","cifUrl":"https://alphafold.ebi.ac.uk/files/other.cif"}]
            """;
        var requested = new List<Uri>();
        using var http = new HttpClient(new RespondingHandler(request =>
        {
            requested.Add(request.RequestUri!);
            return Text(metadata);
        }));
        var candidate = Assert.IsType<CandidateSourceRecord>(
            ExternalSourceExchange.IdentifyExactReference(SourceRouteKind.AlphaFold, "AF-P12345-F1"));
        using var directory = new TemporaryDirectory();

        await Assert.ThrowsAsync<InvalidDataException>(() =>
            new ExternalSourceExchange(http).RetrieveAsync(candidate, directory.Path, CancellationToken.None));
        Assert.Single(requested);
        Assert.Equal("/api/prediction/P12345", requested[0].AbsolutePath);
    }

    private static HttpResponseMessage Bytes(byte[] bytes) => new(HttpStatusCode.OK)
    {
        Content = new ByteArrayContent(bytes)
    };

    private static HttpResponseMessage Text(string value) => new(HttpStatusCode.OK)
    {
        Content = new StringContent(value, Encoding.UTF8, "application/json")
    };

    private sealed class RespondingHandler(Func<HttpRequestMessage, HttpResponseMessage> respond) : HttpMessageHandler
    {
        protected override Task<HttpResponseMessage> SendAsync(
            HttpRequestMessage request, CancellationToken cancellationToken) => Task.FromResult(respond(request));
    }
}

internal sealed class TemporaryDirectory : IDisposable
{
    public string Path { get; } = Directory.CreateTempSubdirectory("pim-slice1-").FullName;
    public void Dispose() => Directory.Delete(Path, recursive: true);
}
