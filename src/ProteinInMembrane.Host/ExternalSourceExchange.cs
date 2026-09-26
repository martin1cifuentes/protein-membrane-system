using System.Collections.Immutable;
using System.Net;
using System.Net.Http.Json;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.RegularExpressions;
using ProteinInMembrane.Host.ProteinInMembraneSystem;

namespace ProteinInMembrane.Host;

public sealed record SourceSearchResult(
    ImmutableArray<CandidateSourceRecord> Candidates,
    ImmutableArray<string> UnavailableRoutes);

// Search hits are candidates. This exchange retrieves original source bytes;
// Protein Preparation decides whether the selected coordinates are usable.
public partial class ExternalSourceExchange
{
    private const long MaximumCoordinateBytes = 100_000_000;
    private const long MaximumPredictionEvidenceBytes = 100_000_000;
    private readonly HttpClient _http;

    public ExternalSourceExchange(HttpClient http)
    {
        _http = http;
        _http.Timeout = TimeSpan.FromSeconds(45);
        // AlphaFold DB rejects anonymous HTTP clients without a User-Agent.
        // Identify this bounded local source retriever on every provider request.
        _http.DefaultRequestHeaders.UserAgent.ParseAdd("ProteinInMembraneSystem/0.1");
    }

    // An exact reference is independent of candidate discovery. Require the
    // AlphaFold model/fragment identity here; an accession may name several.
    public static CandidateSourceRecord? IdentifyExactReference(SourceRouteKind kind, string identifier)
    {
        if (string.IsNullOrWhiteSpace(identifier)) return null;
        identifier = identifier.Trim().ToUpperInvariant();
        return kind switch
        {
            SourceRouteKind.Rcsb when PdbId().IsMatch(identifier) =>
                new CandidateSourceRecord($"rcsb:{identifier}", $"PDB {identifier}", kind,
                    "Researcher-supplied exact RCSB PDB entry reference",
                    ImmutableArray<string>.Empty, null),
            SourceRouteKind.AlphaFold when AlphaFoldModelId().IsMatch(identifier) =>
                new CandidateSourceRecord($"alphafold:{identifier}", $"AlphaFold DB {identifier}", kind,
                    "Researcher-supplied exact AlphaFold DB model or fragment reference",
                    ImmutableArray.Create("Predicted structure; not experimental validation."), null),
            _ => null
        };
    }

    public async Task<SourceSearchResult> SearchAsync(string query, CancellationToken cancellationToken)
    {
        query = query.Trim();
        if (query.Length is < 2 or > 160)
            return new SourceSearchResult(ImmutableArray<CandidateSourceRecord>.Empty,
                ImmutableArray.Create("Enter a descriptive query between 2 and 160 characters."));

        var rcsbTask = SearchRcsbAsync(query, cancellationToken);
        var predictedTask = SearchAlphaFoldAsync(query, cancellationToken);
        await Task.WhenAll(rcsbTask, predictedTask);
        var unavailable = ImmutableArray.CreateBuilder<string>();
        if (rcsbTask.Result.Issue is { } rcsbIssue) unavailable.Add(rcsbIssue);
        if (predictedTask.Result.Issue is { } predictedIssue) unavailable.Add(predictedIssue);
        return new SourceSearchResult(rcsbTask.Result.Candidates.AddRange(predictedTask.Result.Candidates),
            unavailable.ToImmutable());
    }

    public async Task<StructuralSource> RetrieveAsync(
        CandidateSourceRecord candidate,
        string targetDirectory,
        CancellationToken cancellationToken)
    {
        Directory.CreateDirectory(targetDirectory);
        if (candidate.Kind == SourceRouteKind.Rcsb)
        {
            if (!candidate.Id.StartsWith("rcsb:", StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException("The selected RCSB reference is invalid.");
            var accession = candidate.Id["rcsb:".Length..];
            if (!PdbId().IsMatch(accession))
                throw new InvalidDataException("The selected RCSB identifier is invalid.");
            var address = new Uri($"https://files.rcsb.org/download/{accession.ToUpperInvariant()}.cif");
            return await SaveRemoteSourceAsync(candidate, address, targetDirectory, cancellationToken);
        }
        if (candidate.Kind == SourceRouteKind.AlphaFold)
        {
            if (!candidate.Id.StartsWith("alphafold:", StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException("The selected AlphaFold reference is invalid.");
            var selectedId = candidate.Id["alphafold:".Length..];
            var modelMatch = AlphaFoldModelId().Match(selectedId);
            var accession = modelMatch.Success ? modelMatch.Groups[1].Value : selectedId;
            if (!UniProtId().IsMatch(accession))
                throw new InvalidDataException("The selected AlphaFold accession is invalid.");
            var metadataAddress = new Uri($"https://alphafold.ebi.ac.uk/api/prediction/{Uri.EscapeDataString(accession)}");
            using var metadata = await _http.GetAsync(metadataAddress, cancellationToken);
            metadata.EnsureSuccessStatusCode();
            using var document = await JsonDocument.ParseAsync(await metadata.Content.ReadAsStreamAsync(cancellationToken), cancellationToken: cancellationToken);
            if (document.RootElement.ValueKind != JsonValueKind.Array || document.RootElement.GetArrayLength() == 0)
                throw new InvalidDataException("The AlphaFold record has no identified prediction.");
            var predictions = document.RootElement.EnumerateArray()
                .Where(item => item.ValueKind == JsonValueKind.Object)
                .Where(item => !modelMatch.Success ||
                    string.Equals(PredictionRecordId(item), selectedId, StringComparison.OrdinalIgnoreCase))
                .ToArray();
            if (predictions.Length != 1)
                throw new InvalidDataException("The AlphaFold source does not resolve to exactly one identified prediction; select an exact model or fragment.");
            var prediction = predictions[0];
            var recordId = PredictionRecordId(prediction);
            if (string.IsNullOrWhiteSpace(recordId))
                throw new InvalidDataException("The AlphaFold prediction lacks an exact record identity.");
            var coordinateUrl = OptionalPropertyText(prediction, "cifUrl");
            if (!PermittedAlphaFoldAsset(coordinateUrl, out var address))
                throw new InvalidDataException("The AlphaFold record did not supply a permitted mmCIF address.");
            var identifiedCandidate = candidate with { Id = $"alphafold:{recordId}", Label = $"AlphaFold DB {recordId}" };
            var source = await SaveRemoteSourceAsync(identifiedCandidate, address!, targetDirectory, cancellationToken);
            var paeUrl = OptionalPropertyText(prediction, "paeDocUrl");
            string? paePath = null;
            string? paeSha256 = null;
            PaeAcquisitionStanding paeStanding;
            string? paeReason = null;
            if (paeUrl is null)
            {
                paeStanding = PaeAcquisitionStanding.Unavailable;
                paeReason = "The identified AlphaFold prediction did not supply a PAE data address.";
            }
            else if (!PermittedAlphaFoldAsset(paeUrl, out var paeAddress))
            {
                paeStanding = PaeAcquisitionStanding.Unavailable;
                paeReason = "The identified prediction supplied a PAE address outside the permitted source.";
            }
            else
            {
                try
                {
                    (paePath, paeSha256) = await SavePredictionEvidenceAsync(paeAddress!, targetDirectory, cancellationToken);
                    paeStanding = PaeAcquisitionStanding.Available;
                }
                catch (Exception exception) when (!cancellationToken.IsCancellationRequested &&
                    exception is HttpRequestException or TaskCanceledException or IOException or InvalidDataException)
                {
                    paeStanding = PaeAcquisitionStanding.Unavailable;
                    paeReason = $"PAE data for the identified prediction could not be acquired: {exception.Message}";
                }
            }
            var version = OptionalPropertyText(prediction, "latestVersion") ??
                OptionalPropertyText(prediction, "modelVersion");
            return source with
            {
                Prediction = new PredictionEvidenceAsset(recordId, version,
                    OptionalPropertyInt(prediction, "sequenceStart"),
                    OptionalPropertyInt(prediction, "sequenceEnd"),
                    address!.ToString(), source.Sha256, paeUrl, paePath, paeSha256,
                    paeStanding, paeReason),
                Provenance = $"{source.Provenance}; prediction metadata {metadataAddress}"
            };
        }
        throw new InvalidDataException("The selected source kind has no established retrieval route.");
    }

    // OPM's oriented PDB asset is optional reference evidence, not an orientation
    // for the prepared construct. The exact correspondence still belongs to
    // Placement Assessment; missing or uninterpretable records stay unavailable.
    public virtual async Task<OpmReferenceRecord?> TryRetrieveOpmReferenceAsync(
        string pdbAccession,
        string targetDirectory,
        CancellationToken cancellationToken)
    {
        if (!PdbId().IsMatch(pdbAccession)) return null;
        var accession = pdbAccession.ToLowerInvariant();
        var address = new Uri($"https://biomembhub.org/shared/opm-assets/pdb/{accession}.pdb");
        try
        {
            using var response = await _http.GetAsync(address, HttpCompletionOption.ResponseHeadersRead, cancellationToken);
            if (response.StatusCode == HttpStatusCode.NotFound) return null;
            response.EnsureSuccessStatusCode();
            if (response.Content.Headers.ContentLength is > MaximumCoordinateBytes) return null;
            Directory.CreateDirectory(targetDirectory);
            var path = Path.Combine(targetDirectory, $"opm-{accession}-{Guid.NewGuid():N}.pdb");
            await using (var target = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            await using (var input = await response.Content.ReadAsStreamAsync(cancellationToken))
            {
                var buffer = new byte[64 * 1024];
                long total = 0;
                int count;
                while ((count = await input.ReadAsync(buffer.AsMemory(), cancellationToken)) > 0)
                {
                    total += count;
                    if (total > MaximumCoordinateBytes) throw new InvalidDataException("The OPM coordinate record exceeds the accepted local intake size.");
                    await target.WriteAsync(buffer.AsMemory(0, count), cancellationToken);
                }
            }
            var hasCoordinates = File.ReadLines(path).Any(line => line.StartsWith("ATOM  ", StringComparison.Ordinal) ||
                line.StartsWith("HETATM", StringComparison.Ordinal));
            if (!hasCoordinates) return null;
            await using var stream = File.OpenRead(path);
            var hash = Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken)).ToLowerInvariant();
            return new OpmReferenceRecord(pdbAccession.ToUpperInvariant(), address.ToString(), path, hash,
                null, ImmutableArray<ChainSelection>.Empty, ImmutableArray<string>.Empty,
                ImmutableArray<string>.Empty, "OPM oriented PDB; exact assembly and membrane context not established",
                null, null);
        }
        catch (Exception exception) when (exception is HttpRequestException or TaskCanceledException or IOException or InvalidDataException)
        {
            return null;
        }
    }

    private async Task<(ImmutableArray<CandidateSourceRecord> Candidates, string? Issue)> SearchRcsbAsync(string query, CancellationToken cancellationToken)
    {
        try
        {
            var body = new
            {
                query = new { type = "terminal", service = "full_text", parameters = new { value = query } },
                return_type = "entry",
                request_options = new { paginate = new { start = 0, rows = 12 }, results_content_type = new[] { "experimental" } }
            };
            using var response = await _http.PostAsJsonAsync("https://search.rcsb.org/rcsbsearch/v2/query", body, cancellationToken);
            if (response.StatusCode == HttpStatusCode.NoContent)
                return (ImmutableArray<CandidateSourceRecord>.Empty, null);
            response.EnsureSuccessStatusCode();
            using var document = await JsonDocument.ParseAsync(await response.Content.ReadAsStreamAsync(cancellationToken), cancellationToken: cancellationToken);
            if (!document.RootElement.TryGetProperty("result_set", out var hits) || hits.ValueKind != JsonValueKind.Array)
                return (ImmutableArray<CandidateSourceRecord>.Empty, "RCSB returned an uninterpretable candidate account.");
            var builder = ImmutableArray.CreateBuilder<CandidateSourceRecord>();
            foreach (var hit in hits.EnumerateArray())
            {
                if (!hit.TryGetProperty("identifier", out var identifier))
                    continue;
                var id = identifier.GetString();
                if (id is null || !PdbId().IsMatch(id))
                    continue;
                builder.Add(new CandidateSourceRecord($"rcsb:{id.ToUpperInvariant()}", $"PDB {id.ToUpperInvariant()}",
                    SourceRouteKind.Rcsb, "RCSB PDB experimental entry search", ImmutableArray<string>.Empty, null));
            }
            return (builder.ToImmutable(), null);
        }
        catch (Exception exception) when (exception is HttpRequestException or TaskCanceledException or JsonException)
        {
            return (ImmutableArray<CandidateSourceRecord>.Empty, "RCSB candidate discovery is unavailable; this does not indicate a structural defect.");
        }
    }

    private async Task<(ImmutableArray<CandidateSourceRecord> Candidates, string? Issue)> SearchAlphaFoldAsync(string query, CancellationToken cancellationToken)
    {
        try
        {
            var address = $"https://www.ebi.ac.uk/ebisearch/ws/rest/alphafold?query={Uri.EscapeDataString(query)}&format=json&size=12";
            using var response = await _http.GetAsync(address, cancellationToken);
            response.EnsureSuccessStatusCode();
            using var document = await JsonDocument.ParseAsync(await response.Content.ReadAsStreamAsync(cancellationToken), cancellationToken: cancellationToken);
            if (!document.RootElement.TryGetProperty("entries", out var hits) || hits.ValueKind != JsonValueKind.Array)
                return (ImmutableArray<CandidateSourceRecord>.Empty, "EBI Search returned an uninterpretable prediction account.");
            var builder = ImmutableArray.CreateBuilder<CandidateSourceRecord>();
            foreach (var hit in hits.EnumerateArray())
            {
                if (!hit.TryGetProperty("id", out var identifier))
                    continue;
                var raw = identifier.GetString();
                var model = raw is null ? Match.Empty : AlphaFoldModelId().Match(raw);
                var id = model.Success ? model.Value : raw;
                var accession = model.Success ? model.Groups[1].Value : id;
                if (id is null || accession is null || !UniProtId().IsMatch(accession))
                    continue;
                builder.Add(new CandidateSourceRecord($"alphafold:{id}", $"AlphaFold DB {raw}",
                    SourceRouteKind.AlphaFold, "AlphaFold DB prediction search through EMBL-EBI", ImmutableArray.Create("Predicted structure; not experimental validation."), null));
            }
            return (builder.ToImmutable(), null);
        }
        catch (Exception exception) when (exception is HttpRequestException or TaskCanceledException or JsonException)
        {
            return (ImmutableArray<CandidateSourceRecord>.Empty, "AlphaFold DB candidate discovery is unavailable; this does not indicate a structural defect.");
        }
    }

    private async Task<StructuralSource> SaveRemoteSourceAsync(
        CandidateSourceRecord candidate,
        Uri address,
        string targetDirectory,
        CancellationToken cancellationToken)
    {
        using var response = await _http.GetAsync(address, HttpCompletionOption.ResponseHeadersRead, cancellationToken);
        response.EnsureSuccessStatusCode();
        if (response.Content.Headers.ContentLength is > MaximumCoordinateBytes)
            throw new InvalidDataException("The selected coordinate source exceeds the accepted local intake size.");

        var path = Path.Combine(targetDirectory, $"{Guid.NewGuid():N}.cif");
        await using (var target = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None))
        await using (var input = await response.Content.ReadAsStreamAsync(cancellationToken))
        {
            var buffer = new byte[64 * 1024];
            long total = 0;
            int count;
            while ((count = await input.ReadAsync(buffer.AsMemory(), cancellationToken)) > 0)
            {
                total += count;
                if (total > MaximumCoordinateBytes)
                    throw new InvalidDataException("The coordinate transfer exceeded the accepted local intake size.");
                await target.WriteAsync(buffer.AsMemory(0, count), cancellationToken);
            }
        }
        await using var stream = File.OpenRead(path);
        var hash = Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken)).ToLowerInvariant();
        return new StructuralSource(candidate.Id, candidate.Kind, $"{candidate.Provenance}; {address}",
            path, hash, candidate.Id.Split(':').Last(), candidate.Label);
    }

    private async Task<(string Path, string Sha256)> SavePredictionEvidenceAsync(
        Uri address, string targetDirectory, CancellationToken cancellationToken)
    {
        using var response = await _http.GetAsync(address, HttpCompletionOption.ResponseHeadersRead, cancellationToken);
        response.EnsureSuccessStatusCode();
        if (response.Content.Headers.ContentLength is > MaximumPredictionEvidenceBytes)
            throw new InvalidDataException("The identified PAE data exceeds the bounded local evidence size.");
        Directory.CreateDirectory(targetDirectory);
        var path = Path.Combine(targetDirectory, $"pae-{Guid.NewGuid():N}.json");
        await using (var target = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None))
        await using (var input = await response.Content.ReadAsStreamAsync(cancellationToken))
        {
            var buffer = new byte[64 * 1024];
            long total = 0;
            int count;
            while ((count = await input.ReadAsync(buffer.AsMemory(), cancellationToken)) > 0)
            {
                total += count;
                if (total > MaximumPredictionEvidenceBytes)
                    throw new InvalidDataException("The PAE transfer exceeded its bounded local evidence size.");
                await target.WriteAsync(buffer.AsMemory(0, count), cancellationToken);
            }
        }
        await using var stream = File.OpenRead(path);
        return (path, Convert.ToHexString(await SHA256.HashDataAsync(stream, cancellationToken)).ToLowerInvariant());
    }

    private static string? PredictionRecordId(JsonElement prediction)
        => OptionalPropertyText(prediction, "modelEntityId") ?? OptionalPropertyText(prediction, "entryId");

    private static string? OptionalPropertyText(JsonElement item, string name)
    {
        if (!item.TryGetProperty(name, out var value)) return null;
        return value.ValueKind switch
        {
            JsonValueKind.String => value.GetString(),
            JsonValueKind.Number => value.GetRawText(),
            _ => null
        };
    }

    private static int? OptionalPropertyInt(JsonElement item, string name)
        => item.TryGetProperty(name, out var value) && value.ValueKind == JsonValueKind.Number &&
            value.TryGetInt32(out var result) ? result : null;

    private static bool PermittedAlphaFoldAsset(string? value, out Uri? address)
    {
        address = null;
        return Uri.TryCreate(value, UriKind.Absolute, out address) &&
            address.Scheme == Uri.UriSchemeHttps &&
            address.Host is "alphafold.ebi.ac.uk" or "alphafold.com";
    }

    [GeneratedRegex("^[0-9A-Za-z]{4}$", RegexOptions.CultureInvariant)]
    private static partial Regex PdbId();

    [GeneratedRegex("^[A-Za-z0-9]{6,10}$", RegexOptions.CultureInvariant)]
    private static partial Regex UniProtId();

    [GeneratedRegex("^AF-([A-Za-z0-9]{6,10})-F[1-9][0-9]*$", RegexOptions.CultureInvariant | RegexOptions.IgnoreCase)]
    private static partial Regex AlphaFoldModelId();
}
