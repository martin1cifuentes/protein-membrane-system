using System.Net;
using Microsoft.Extensions.FileProviders;
using ProteinInMembrane.Host;
using ProductRoot = ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem;

var originText = Environment.GetEnvironmentVariable("ASPNETCORE_URLS") ?? "http://127.0.0.1:4185";
if (!Uri.TryCreate(originText, UriKind.Absolute, out var origin) ||
    origin.Scheme != Uri.UriSchemeHttp || origin.Host != IPAddress.Loopback.ToString() ||
    origin.Port is < 1024 or > 65535 || origin.AbsolutePath != "/" ||
    origin.Query.Length != 0 || origin.Fragment.Length != 0)
    throw new InvalidOperationException("The local host must use one HTTP address on 127.0.0.1 and an unprivileged port.");

var python = RequiredEnvironmentPath("PIM_WORKER_PYTHON");
var ownerSourceRoot = RequiredEnvironmentPath("PIM_OWNER_SOURCE_ROOT");
var workspaceRoot = RequiredEnvironmentPath("PIM_WORKSPACE_ROOT");
var browserRoot = Path.Combine(AppContext.BaseDirectory, "wwwroot");
if (!File.Exists(Path.Combine(browserRoot, "index.html")))
    throw new InvalidOperationException("The local browser bundle is absent. Build the product before starting it.");

var builder = WebApplication.CreateBuilder(new WebApplicationOptions
{
    Args = args,
    WebRootPath = browserRoot
});
builder.WebHost.ConfigureKestrel(options =>
{
    options.Listen(IPAddress.Loopback, origin.Port);
    options.Limits.MaxRequestBodySize = 101_000_000;
});

// Source URLs are selected by this product. Do not follow a provider redirect
// to an unreviewed host or local address as part of structure retrieval.
using var remoteClient = new HttpClient(new HttpClientHandler { AllowAutoRedirect = false });
var sources = new ExternalSourceExchange(remoteClient);
var worker = new ScientificWorkerExchange(python, ownerSourceRoot);
var system = new ProductRoot(worker, sources, workspaceRoot,
    Environment.GetEnvironmentVariable("PIM_POLICY_CATALOGUE") ?? string.Empty,
    Environment.GetEnvironmentVariable("PIM_PPM_EXECUTABLE") ?? string.Empty,
    Environment.GetEnvironmentVariable("PIM_PACKMOL_EXECUTABLE") ?? string.Empty);

var app = builder.Build();
app.UseExceptionHandler(error => error.Run(async context =>
{
    context.Response.StatusCode = StatusCodes.Status500InternalServerError;
    await context.Response.WriteAsJsonAsync(new { reason = "The local product encountered an unexpected error." });
}));
app.Use(async (context, next) =>
{
    if (!context.Request.Host.Host.Equals(IPAddress.Loopback.ToString(), StringComparison.OrdinalIgnoreCase) ||
        context.Request.Host.Port != origin.Port ||
        !IPAddress.IsLoopback(context.Connection.RemoteIpAddress ?? IPAddress.None))
    {
        context.Response.StatusCode = StatusCodes.Status403Forbidden;
        return;
    }
    context.Response.Headers["X-Content-Type-Options"] = "nosniff";
    context.Response.Headers["X-Frame-Options"] = "DENY";
    context.Response.Headers["Referrer-Policy"] = "no-referrer";
    context.Response.Headers.CacheControl = "no-store";
    await next();
});

BrowserExchange.Map(app, system, origin);
var browserFiles = new PhysicalFileProvider(browserRoot);
app.UseDefaultFiles(new DefaultFilesOptions { FileProvider = browserFiles });
app.UseStaticFiles(new StaticFileOptions { FileProvider = browserFiles });
app.MapFallback((HttpContext context) =>
    context.Request.Path.StartsWithSegments("/api", StringComparison.OrdinalIgnoreCase)
        ? Results.NotFound()
        : Results.File(Path.Combine(browserRoot, "index.html"), "text/html"));

await app.RunAsync();

static string RequiredEnvironmentPath(string name)
{
    var path = Environment.GetEnvironmentVariable(name);
    if (string.IsNullOrWhiteSpace(path))
        throw new InvalidOperationException($"The local prerequisite {name} is not configured.");
    return Path.GetFullPath(path);
}
