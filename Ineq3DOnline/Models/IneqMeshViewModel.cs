using MeshData;
using Newtonsoft.Json.Serialization;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshViewModel
    {
        private static string[] sampleUFuncs = {
            "sin(6 * x - t) + sin(6 * y - t) + sin(6 * z - t)",
            "x * y * (z - sin(t)^2)",
            "sin(6 * x - 0.1*t^2) + sin(6 * y - t) + sin(6 * z - t)",
            "(x-sin(t))^2+y^2+z^2",
            "(x-sin(t))^2+cos(t)*y^2+z^2",
            "x * y ^2 * (z - sin(t)) +0.4*y*(x-cos(2*t))^3+z",
            "sin(5*((x-sin(t))^2+cos(t)*y^2+z^2))",
            "sin(4*x-4*sin(t))+sin(4*y)+sin(4*z)",
            "sin(t)*(x^2+y^2+z^2)+(1-sin(t))*(x*y*z)",
            "sin(t)*(z-x^2+y^2)+(1-sin(t))*(x-y^2-z^2)"
        } ;

        public int SampleUFuncIndex { get; set; }

        public static IneqMeshViewModel DefaultModel(int sampleU = -1)
        {
            if (sampleU == -1)
            {
                sampleU = new Random().Next(sampleUFuncs.Length);
            }
            else
            {
                sampleU = sampleU % sampleUFuncs.Length;
            }

            return new IneqMeshViewModel
            {
                Name = "",
                Formula = "x^2 + y^2 + z^2 < 1",
                SampleUFuncIndex = sampleU,
                MaxDivisionCount = 12,
                X0 = -1,
                Y0 = -1,
                Z0 = -1,
                X1 = 1,
                Y1 = 1,
                Z1 = 1,
                UFunc = sampleUFuncs[sampleU]
            };
        }

        [Required(ErrorMessage = "The Inequalities are required.")]
        public string Formula { get; set; }

        public string Name { get; set; }

        [Required(ErrorMessage = "The Division count is required.")]
        [Range(5, 25, ErrorMessage = "Value for {0} must be between {1} and {2}.")]
        [Display(Name = "Division count")]
        public int MaxDivisionCount { get; set; }

        public double X0 { get; set; }
        public double Y0 { get; set; }
        public double Z0 { get; set; }
        public double X1 { get; set; }
        public double Y1 { get; set; }
        public double Z1 { get; set; }

        [Display(Name = "Improve mesh quality")]
        public bool Quality { get; set; }

        [Display(Name = "Improve approximation by curvature")]
        public bool CurvatureQuality { get; set; }

        //[Required(ErrorMessage = "The Function is required.")]
        [Display(Name = "u(x,y,z,t) = ")]
        public string UFunc { get; set; }

        public IneqMesh IneqMesh { get; set; }

        public string PLY { get; set; }

        public void SetIneqMesh()
        {
            if (Formula.StartsWith("[advancedX]") || (HttpContext.Current.Request.IsLocal && Formula.StartsWith("[advanced]")))
            {
                IneqMesh = DynamicCodeExecutor.Execute(Formula.Replace("[advancedX]", "").Replace("[advanced]", ""));
            }
            else
            {
                IneqMesh = new IneqMesh
                {
                    X0 = X0,
                    Y0 = Y0,
                    Z0 = Z0,
                    X1 = X1,
                    Y1 = Y1,
                    Z1 = Z1,
                    D = Math.Max(Math.Max(X1 - X0, Y1 - Y0), Z1 - Z0) / MaxDivisionCount,
                    Boxed = true,
                    IneqTree = IneqTreeParser.FromFormula(Formula)
                };
            }
        }

        public void CreateMesh(bool createPLY = true) 
        {
            if (IneqMesh == null)
            {
                return;
            }

            IneqMesh.Create();

            if (Quality)
            {
                IneqMeshTools.CheckQuality(IneqMesh);
                IneqMeshTools.CheckBoundaryQuality(IneqMesh);
            }

            if (CurvatureQuality)
            {
                IneqMeshTools.CheckCurvatureQuality(IneqMesh);
            }
            
            IneqMesh.DeleteLonelyPoints();

            if (createPLY)
            {
                PLY = PLYTools.GetPLY(IneqMesh);
            }
        }

        public void Save(string name, bool saveJSON, bool savePLY)
        {
            Name = name;

            if (saveJSON)
            {

                string path = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/Samples"), name + ".json");
                System.IO.File.WriteAllText(path, JsonConvert.SerializeObject(new
                {
                    Name,
                    X0,
                    Y0,
                    Z0,
                    X1,
                    Y1,
                    Z1,
                    MaxDivisionCount,
                    Quality,
                    CurvatureQuality,
                    Formula = Formula.Replace("[advancedX]", "[advanced]")
                }, Formatting.None, new JsonSerializerSettings
                {
                    ContractResolver = new CamelCasePropertyNamesContractResolver(),
                    DateFormatHandling = DateFormatHandling.IsoDateFormat,
                    NullValueHandling = NullValueHandling.Ignore
                }));
            }

            if (savePLY &&!String.IsNullOrEmpty(PLY))
            {
                string path = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/Samples"), name + ".ply");
                System.IO.File.WriteAllText(path, PLY);
            }
        }

        public IEnumerable<string> DataSamples
        {
            get
            {
                string path = HttpContext.Current.Server.MapPath("~/Samples");
                return System.IO.Directory.GetFiles(path, "*.json").Select(fn => System.IO.Path.GetFileNameWithoutExtension(fn));
            }
        }
    }
}