using MeshData;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshViewModel
    {

        private static string[] sampleFormulas =
        {   @"z < (sin(6*x)+sin(6*y))/4",
            @"(
	x^2 + y^2 + z^2 < 0.75^2
) ||
(
	x^2 + y^2 + (z-0.75)^2 < 0.25^2
)",
            @"(
    x^2 + y^2 + z^2 < 1.25 ^ 2
) &&
x^2+y^2 > 0.5^2 &&
x^2+z^2 > 0.5^2 &&
z^2+y^2 > 0.5^2",

            @"(
    x^2 + y^2 + z^2 < 1 &&
    x + y < z^2
) ||
(
    x^2 + z^2 < 0.25
)",
            @"(
	x^2+y^2+z^2<1 &&
	z<0+(cos(5*x)+cos(5*y))/3
) ||
x^2+y^2+(z-.8)^2<0.25"
        };


        public static IneqMeshViewModel DefaultModel()
        {
            int sample = new Random().Next(sampleFormulas.Length);

            return new IneqMeshViewModel
            {
                Formula = sampleFormulas[sample],
                MaxDivisionCount = 12,
                X0 = -1,
                Y0 = -1,
                Z0 = -1,
                X1 = 1,
                Y1 = 1,
                Z1 = 1
            };
        }

        [Required(ErrorMessage = "The Inequalities are required.")]
        public string Formula { get; set; }

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
        
        public IneqMesh IneqMesh { get; set; }

        public void SetIneqMesh()
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

        public string PLY { get; set; }
    }
}