using System.Buffers.Binary;
using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text;

namespace ProteinInMembrane.Host.ProteinInMembraneSystem;

/// <summary>
/// Identifies the exact declared optional procedure, independently of JSON property order,
/// serializer options, or the evidence record that qualifies the procedure.
/// </summary>
public static class EquilibrationProtocolFingerprint
{
    public static bool TryCompute(EquilibrationProtocol? protocol, out string digest)
    {
        digest = string.Empty;
        if (protocol is null) return false;
        try
        {
            digest = Compute(protocol);
            return true;
        }
        catch (Exception exception) when (exception is ArgumentException or NullReferenceException or
                                     InvalidOperationException)
        {
            return false;
        }
    }

    public static string Compute(EquilibrationProtocol protocol)
    {
        ArgumentNullException.ThrowIfNull(protocol);
        using var bytes = new MemoryStream();
        PutString(bytes, "equilibration-protocol-v1");
        PutString(bytes, protocol.Id);
        PutDouble(bytes, protocol.TargetTemperatureKelvin);
        PutInt(bytes, protocol.RandomSeed);
        PutArray(bytes, protocol.Stages, PutStage);
        PutStage(bytes, protocol.ExtensionWindow);
        PutInt(bytes, protocol.MaximumExtensions);
        PutInt(bytes, protocol.MaximumSampleCount);
        PutArray(bytes, protocol.RequiredObservations, PutString);
        PutArray(bytes, protocol.Observables, PutObservable);
        PutArray(bytes, protocol.SufficiencyRules, PutRule);
        PutString(bytes, protocol.ComparisonBasis);
        return Convert.ToHexString(SHA256.HashData(bytes.ToArray())).ToLowerInvariant();
    }

    private static void PutStage(Stream target, EquilibrationStageControl stage)
    {
        ArgumentNullException.ThrowIfNull(stage);
        PutString(target, stage.Name);
        PutInt(target, stage.Steps);
        PutDouble(target, stage.TimestepPicoseconds);
        PutDouble(target, stage.TemperatureKelvin);
        PutNullableDouble(target, stage.PressureBar);
        PutString(target, stage.PressureMode);
        PutNullableInt(target, stage.BarostatFrequencySteps);
        PutNullableDouble(target, stage.SurfaceTensionBarNm);
        PutDouble(target, stage.FrictionPerPicosecond);
        PutDouble(target, stage.ProteinRestraintKjMolNm2);
        PutDouble(target, stage.LipidRestraintKjMolNm2);
        PutInt(target, stage.ReportIntervalSteps);
    }

    private static void PutObservable(Stream target, EquilibrationObservable item)
    {
        PutString(target, item.Name);
        PutString(target, item.Unit);
        PutString(target, item.Scope);
        PutString(target, item.Category);
        PutString(target, item.Method);
        PutString(target, item.AtomSelector);
        PutString(target, item.ComparisonSelector);
        PutString(target, item.ThirdSelector);
        PutNullableDouble(target, item.DistanceCutoffAngstrom);
    }

    private static void PutRule(Stream target, EquilibrationSufficiencyRule rule)
    {
        PutString(target, rule.ObservableName);
        PutInt(target, rule.BlockSizeSamples);
        PutInt(target, rule.MinimumEffectiveBlocks);
        PutDouble(target, rule.MaximumAbsoluteFirstVsLastBlockMeanDifference);
        PutDouble(target, rule.MaximumAbsoluteLagOneBlockCorrelation);
    }

    private static void PutArray<T>(Stream target, ImmutableArray<T> values, Action<Stream, T> put)
    {
        if (values.IsDefault) throw new ArgumentException("A protocol array is absent.");
        PutInt(target, values.Length);
        foreach (var value in values) put(target, value);
    }

    private static void PutString(Stream target, string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        var utf8 = Encoding.UTF8.GetBytes(value);
        PutInt(target, utf8.Length);
        target.Write(utf8);
    }

    private static void PutNullableInt(Stream target, int? value)
    {
        target.WriteByte(value.HasValue ? (byte)1 : (byte)0);
        if (value.HasValue) PutInt(target, value.Value);
    }

    private static void PutNullableDouble(Stream target, double? value)
    {
        target.WriteByte(value.HasValue ? (byte)1 : (byte)0);
        if (value.HasValue) PutDouble(target, value.Value);
    }

    private static void PutInt(Stream target, int value)
    {
        Span<byte> buffer = stackalloc byte[sizeof(int)];
        BinaryPrimitives.WriteInt32LittleEndian(buffer, value);
        target.Write(buffer);
    }

    private static void PutDouble(Stream target, double value)
    {
        Span<byte> buffer = stackalloc byte[sizeof(long)];
        BinaryPrimitives.WriteInt64LittleEndian(buffer, BitConverter.DoubleToInt64Bits(value));
        target.Write(buffer);
    }
}
