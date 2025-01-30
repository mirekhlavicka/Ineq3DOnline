using MeshData;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;
using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Web.SessionState;

namespace Ineq3DOnline
{
    public class DynamicCodeExecutor
    {
        private static string codeTemplate = @"
            using System;
            using System.Collections.Generic;
            using MeshData;
            using Ineq3DOnline;

            namespace DynamicNamespace
            {
                public class DynamicClass
                {
                    {0}
                }
            }";

        public static IneqMesh Execute(string code)
        {
            // Create a syntax tree
            var syntaxTree = CSharpSyntaxTree.ParseText(codeTemplate.Replace("{0}", code));

            // Reference necessary assemblies
            var references = AppDomain.CurrentDomain.GetAssemblies()
                .Where(a => !a.IsDynamic && !string.IsNullOrEmpty(a.Location))
                .Select(a => MetadataReference.CreateFromFile(a.Location))
                .ToList();

            // Compile the code
            var compilation = CSharpCompilation.Create("DynamicAssembly")
                .WithOptions(new CSharpCompilationOptions(OutputKind.DynamicallyLinkedLibrary))
                .AddReferences(references)
                .AddSyntaxTrees(syntaxTree);

            using (var ms = new MemoryStream())
            {
                var result = compilation.Emit(ms);
                if (!result.Success)
                {
                    // Handle compilation errors
                    var errors = string.Join(/*Environment.NewLine*/"<br/>", result.Diagnostics
                        .Where(diag => diag.Severity == DiagnosticSeverity.Error)
                        .Select(diag => AdjustLineNumbers(diag.ToString())));
                    throw new Exception($"Compilation failed: {errors}");
                }

                ms.Seek(0, SeekOrigin.Begin);

                // Load the compiled assembly
                var assembly = Assembly.Load(ms.ToArray());

                // Use reflection to invoke the code
                var type = assembly.GetType("DynamicNamespace.DynamicClass");
                var instance = Activator.CreateInstance(type);
                var method = type.GetMethod("GetIneqMesh");

                var resultMessage = method.Invoke(instance, new object[] { });
                return (IneqMesh)resultMessage;
            }
        }

        static string AdjustLineNumbers(string input)
        {
            return Regex.Replace(input, @"\((-?\d+),\s*(-?\d+)\)", match =>
            {
                int n1 = int.Parse(match.Groups[1].Value) - 10;
                int n2 = int.Parse(match.Groups[2].Value);
                return $"({n1}, {n2})";
            });
        }
    }
}