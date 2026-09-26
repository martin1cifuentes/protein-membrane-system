using System.Collections.Immutable;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;
using Microsoft.CodeAnalysis.CSharp.Syntax;

if (args.Length != 1 || !Directory.Exists(args[0]))
{
    Console.Error.WriteLine("Usage: ArchitectureCheck <Host source directory>");
    return 2;
}

var sourceRoot = Path.GetFullPath(args[0]);
var files = Directory.GetFiles(sourceRoot, "*.cs", SearchOption.AllDirectories)
    .Where(path => !Path.GetRelativePath(sourceRoot, path)
        .Split([Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar])
        .Any(segment => segment is "bin" or "obj"))
    .OrderBy(path => path, StringComparer.Ordinal)
    .ToArray();
if (files.Length == 0)
{
    Console.Error.WriteLine("ARCH000: No Host C# source was found.");
    return 1;
}

// Each governed source file has one assigned owner. SourceSearchResult is the
// one passive information type co-located with its source exchange; it is not
// an exchange endpoint. This map is intentionally exact, so a newly added C#
// file must first receive a Source Topology assignment.
const string productNamespace = "ProteinInMembrane.Host.ProteinInMembraneSystem";
var assignments = new Dictionary<string, (string Role, string Namespace)>(StringComparer.Ordinal)
{
    ["Program.cs"] = ("composition", ""),
    ["BrowserExchange.cs"] = ("browserExchange", "ProteinInMembrane.Host"),
    ["ExternalSourceExchange.cs"] = ("sourceExchange", "ProteinInMembrane.Host"),
    ["ScientificWorkerExchange.cs"] = ("workerExchange", "ProteinInMembrane.Host"),
    ["ProteinInMembraneSystem/ProteinInMembraneSystem.cs"] = ("root", productNamespace),
    ["ProteinInMembraneSystem/ProteinInMembraneSystem.Contracts.cs"] = ("shared", productNamespace),
    // Exact immutable protocol identity is shared information, like
    // PreparationPolicyFingerprint in the contracts companion. It grants no
    // stage execution or qualification authority to its callers.
    ["ProteinInMembraneSystem/EquilibrationProtocolFingerprint.cs"] = ("shared", productNamespace),
    ["ProteinInMembraneSystem/ProteinPreparation/ProteinPreparation.cs"] =
        ("proteinPreparation", productNamespace + ".ProteinPreparation"),
    ["ProteinInMembraneSystem/MembraneModelAssessment/MembraneModelAssessment.cs"] =
        ("membraneModelAssessment", productNamespace + ".MembraneModelAssessment"),
    ["ProteinInMembraneSystem/PlacementAssessment/PlacementAssessment.cs"] =
        ("placementAssessment", productNamespace + ".PlacementAssessment"),
    ["ProteinInMembraneSystem/ExplicitPreparation/ExplicitPreparation.cs"] =
        ("explicitPreparation", productNamespace + ".ExplicitPreparation"),
    ["ProteinInMembraneSystem/ExplicitPreparation/Minimization.cs"] =
        ("minimization", productNamespace + ".ExplicitPreparation"),
    ["ProteinInMembraneSystem/ExplicitPreparation/OptionalEquilibrationProcedure.cs"] =
        ("optionalEquilibration", productNamespace + ".ExplicitPreparation"),
    ["ProteinInMembraneSystem/PreparationAssessment/PreparationAssessment.cs"] =
        ("preparationAssessment", productNamespace + ".PreparationAssessment"),
    ["ProteinInMembraneSystem/ConnectedStructuralInspection/ConnectedStructuralInspection.cs"] =
        ("connectedStructuralInspection", productNamespace + ".ConnectedStructuralInspection"),
    ["ProteinInMembraneSystem/CompletedStageExport/CompletedStageExport.cs"] =
        ("completedStageExport", productNamespace + ".CompletedStageExport"),
    ["ProteinInMembraneSystem/LocalRunWorkspace/LocalRunWorkspace.cs"] =
        ("localRunWorkspace", productNamespace + ".LocalRunWorkspace")
};

string RelativeSourcePath(string path) =>
    Path.GetRelativePath(sourceRoot, Path.GetFullPath(path)).Replace('\\', '/');

var trees = files.Select(path => CSharpSyntaxTree.ParseText(File.ReadAllText(path),
    new CSharpParseOptions(LanguageVersion.Preview), path)).ToArray();
var platform = AppContext.GetData("TRUSTED_PLATFORM_ASSEMBLIES") as string;
if (string.IsNullOrWhiteSpace(platform))
{
    Console.Error.WriteLine("ARCH000: The .NET platform reference set is unavailable.");
    return 1;
}
var references = platform.Split(Path.PathSeparator, StringSplitOptions.RemoveEmptyEntries)
    .Select(path => MetadataReference.CreateFromFile(path));
var compilation = CSharpCompilation.Create("ProteinInMembrane.ArchitectureInput", trees, references,
    new CSharpCompilationOptions(OutputKind.ConsoleApplication));

var errors = new HashSet<string>(StringComparer.Ordinal);
var declared = new HashSet<string>(StringComparer.Ordinal);
var present = files.Select(RelativeSourcePath).ToHashSet(StringComparer.Ordinal);
foreach (var path in present.Where(path => !assignments.ContainsKey(path)))
    errors.Add($"ARCH000: C# source has no assigned source role: {path}");
foreach (var path in assignments.Keys.Where(path => !present.Contains(path)))
    errors.Add($"ARCH000: Assigned C# source is absent: {path}");
foreach (var tree in trees)
{
    foreach (var diagnostic in tree.GetDiagnostics().Where(item => item.Severity == DiagnosticSeverity.Error))
        errors.Add($"{diagnostic.Location.GetLineSpan()}: ARCH000: {diagnostic.GetMessage()}");

    var model = compilation.GetSemanticModel(tree, ignoreAccessibility: true);
    var root = tree.GetRoot();
    foreach (var declaration in root.DescendantNodes().OfType<BaseTypeDeclarationSyntax>())
        if (model.GetDeclaredSymbol(declaration) is INamedTypeSymbol type)
        {
            declared.Add(type.ToDisplayString());
            var source = RelativeSourcePath(tree.FilePath);
            if (!assignments.TryGetValue(source, out var assignment) ||
                type.ContainingNamespace.ToDisplayString() != assignment.Namespace ||
                RoleOf(type) == "unassigned")
                errors.Add($"{tree.FilePath}:{tree.GetLineSpan(declaration.Span).StartLinePosition.Line + 1}: " +
                    $"ARCH000: Type does not match its assigned source role and namespace: {type.ToDisplayString()}");
        }

    foreach (var node in root.DescendantNodes().Where(item => item is SimpleNameSyntax or
                 InvocationExpressionSyntax or ObjectCreationExpressionSyntax or
                 ImplicitObjectCreationExpressionSyntax or ElementAccessExpressionSyntax))
    {
        var enclosing = model.GetEnclosingSymbol(node.SpanStart);
        var callerType = enclosing as INamedTypeSymbol ?? enclosing?.ContainingType;
        if (callerType is null) continue; // e.g. namespace imports, not behavioral use.
        // Top-level statements have a compiler-generated containing type.
        // Their assigned source, rather than that synthetic type, owns the call.
        var caller = RelativeSourcePath(tree.FilePath) == "Program.cs" &&
            node.Ancestors().OfType<GlobalStatementSyntax>().Any()
            ? "composition" : RoleOf(callerType);
        var symbols = model.GetSymbolInfo(node);
        ImmutableArray<ISymbol> targets = symbols.Symbol is null
            ? symbols.CandidateSymbols : [symbols.Symbol];
        foreach (var symbol in targets)
        {
            var resolved = symbol is IAliasSymbol alias ? alias.Target : symbol;
            var targetType = resolved switch
            {
                INamedTypeSymbol named => named,
                IMethodSymbol method => method.ContainingType,
                IPropertySymbol property => property.ContainingType,
                IEventSymbol eventSymbol => eventSymbol.ContainingType,
                IFieldSymbol field => field.ContainingType,
                _ => null
            };
            if (targetType is null) continue;
            // Anonymous types and tuples are compiler-generated data shapes,
            // not separately owned behavioral endpoints.
            if (targetType.IsAnonymousType || targetType.IsTupleType) continue;
            if (!AllowedWorkerCapability(caller, targetType))
            {
                var workerPosition = tree.GetLineSpan(node.Span).StartLinePosition;
                errors.Add($"{tree.FilePath}:{workerPosition.Line + 1}:{workerPosition.Character + 1}: " +
                    $"ARCH002: {caller} may not use worker capability " +
                    $"{resolved.ToDisplayString(SymbolDisplayFormat.CSharpErrorMessageFormat)}");
                continue;
            }
            var target = RoleOf(targetType);
            if (Allowed(caller, target)) continue;
            var position = tree.GetLineSpan(node.Span).StartLinePosition;
            errors.Add($"{tree.FilePath}:{position.Line + 1}:{position.Character + 1}: " +
                $"ARCH001: {caller} may not use {target} behavioral endpoint " +
                $"{resolved.ToDisplayString(SymbolDisplayFormat.CSharpErrorMessageFormat)}");
        }
    }
}

var required = new[]
{
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinInMembraneSystem",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ProteinPreparation.ProteinPreparation",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.MembraneModelAssessment.MembraneModelAssessment",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.PlacementAssessment.PlacementAssessment",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.ExplicitPreparation",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.Minimization",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ExplicitPreparation.OptionalEquilibrationProcedure",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.PreparationAssessment.PreparationAssessment",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.ConnectedStructuralInspection.ConnectedStructuralInspection",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.CompletedStageExport.CompletedStageExport",
    "ProteinInMembrane.Host.ProteinInMembraneSystem.LocalRunWorkspace.LocalRunWorkspace"
};
foreach (var type in required.Where(type => !declared.Contains(type)))
    errors.Add($"ARCH000: Required behavioral boundary is absent: {type}");

foreach (var error in errors.OrderBy(item => item, StringComparer.Ordinal))
    Console.Error.WriteLine(error);
if (errors.Count > 0) return 1;
Console.WriteLine($"Architecture boundary check passed for {files.Length} Host source files.");
return 0;

// Ownership is established from the resolved symbol's exact declaration
// origin, not from a namespace or the caller's file alone. In particular, a
// child cannot acquire root authority by changing its namespace declaration.
string RoleOf(INamedTypeSymbol type)
{
    var origins = type.DeclaringSyntaxReferences
        .Select(reference => RelativeSourcePath(reference.SyntaxTree.FilePath))
        .Distinct(StringComparer.Ordinal)
        .ToArray();
    // Roslyn synthesizes the top-level Program type. Its declaration can have
    // a Program.cs syntax reference, so recognize that assigned composition
    // source before applying the ordinary declared-type rule.
    if (type.Name == "Program" && type.ContainingNamespace.IsGlobalNamespace &&
        SymbolEqualityComparer.Default.Equals(type.ContainingAssembly, compilation.Assembly) &&
        present.Contains("Program.cs") && origins.All(origin => origin == "Program.cs"))
        return "composition";
    if (origins.Length > 1) return "unassigned"; // no cross-file partial ownership
    if (origins.Length == 1)
    {
        if (!assignments.TryGetValue(origins[0], out var assignment) ||
            type.ContainingNamespace.ToDisplayString() != assignment.Namespace)
            return "unassigned";
        return origins[0] == "ExternalSourceExchange.cs" &&
            type.ToDisplayString() == "ProteinInMembrane.Host.SourceSearchResult"
            ? "shared" : assignment.Role;
    }
    if (type.ContainingType is not null &&
        SymbolEqualityComparer.Default.Equals(type.ContainingAssembly, compilation.Assembly))
        return RoleOf(type.ContainingType);
    if (SymbolEqualityComparer.Default.Equals(type.ContainingAssembly, compilation.Assembly))
    {
        // A compiler-generated source type without a known owning parent is
        // not allowed to acquire external/shared authority by default.
        return "unassigned";
    }
    return "external";
}

static bool Allowed(string caller, string target)
{
    if (target is "external" or "shared") return true;
    if (target == "sourceExchange") return caller is "root" or "composition" or "sourceExchange";
    if (target == "workerExchange") return caller is "root" or "composition" or "workerExchange";
    if (target == "browserExchange") return caller is "composition" or "browserExchange";
    if (caller == target) return true;
    if (caller == "root" && target is ("proteinPreparation" or "membraneModelAssessment" or
            "placementAssessment" or "explicitPreparation" or "preparationAssessment" or
            "connectedStructuralInspection" or "completedStageExport" or "localRunWorkspace")) return true;
    if (caller == "explicitPreparation" && target is ("minimization" or "optionalEquilibration")) return true;
    return caller is ("composition" or "browserExchange") && target == "root";
}

// Shared information remains broadly readable, but the owner-specific worker
// operations declared in that shared contract are not interchangeable.
static bool AllowedWorkerCapability(string caller, INamedTypeSymbol targetType)
{
    if (targetType.ContainingNamespace.ToDisplayString() !=
        "ProteinInMembrane.Host.ProteinInMembraneSystem") return true;
    // The contract declarations themselves may compose these interfaces;
    // behavioral owners remain restricted to their assigned capability.
    if (caller is "shared" or "root" or "composition" or "workerExchange") return true;
    return targetType.Name switch
    {
        "IScientificWorkerExchange" => false,
        "IProteinPreparationWork" => caller == "proteinPreparation",
        "IMembraneModelAssessmentWork" => caller == "membraneModelAssessment",
        "IPlacementAssessmentWork" => caller == "placementAssessment",
        "IExplicitConstructionWork" => caller == "explicitPreparation",
        "IMinimizationWork" => caller is "explicitPreparation" or "minimization",
        "IOptionalEquilibrationWork" => caller is "explicitPreparation" or "optionalEquilibration",
        "ICompletedStageExportWork" => caller == "completedStageExport",
        _ => true
    };
}
