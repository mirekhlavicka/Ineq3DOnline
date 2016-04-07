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
        {
@"z < (sin(6*x)+sin(6*y))/4",
//------------------------------------------------------------
@"(
	x^2 + y^2 + z^2 < 0.75^2
) ||
(
	x^2 + y^2 + (z-0.75)^2 < 0.25^2
)",
//------------------------------------------------------------
@"(
    x^2 + y^2 + z^2 < 1.25 ^ 2
) &&
x^2+y^2 > 0.5^2 &&
x^2+z^2 > 0.5^2 &&
z^2+y^2 > 0.5^2",
//------------------------------------------------------------
@"(
    x^2 + y^2 + z^2 < 1 &&
    x + y < z^2
) ||
(
    x^2 + z^2 < 0.25
)",
//------------------------------------------------------------
@"(
	x^2+y^2+z^2<1 &&
	z<0+(cos(5*x)+cos(5*y))/3
) ||
x^2+y^2+(z-.8)^2<0.25",
//------------------------------------------------------------
@"abs(x)<0.25 ||
0.5*x^2+y^2+z^2<0.5",
//------------------------------------------------------------
@"8*x^2+8*y^2+z^2<1 ||
8*x^2+y^2+8*z^2<1 ||
x^2+8*y^2+8*z^2<1",
//------------------------------------------------------------

@"16*x^2+16*y^2+z^2<1 ||
90*x^6+y^2+90*z^6<1 ||
x^2+8*y^2+8*z^2<1"

        };

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


        public int SampleIndex { get; set; }
        public int SampleUFuncIndex { get; set; }

        public static IneqMeshViewModel DefaultModel(int sample = -1, int sampleU = -1)
        {
            if (sample == -1)
            {
                sample = new Random().Next(sampleFormulas.Length);
            }
            else
            {
                sample = sample % sampleFormulas.Length;
            }

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
                Formula = sampleFormulas[sample],
                SampleIndex = sample,
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