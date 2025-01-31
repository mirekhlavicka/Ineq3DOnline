﻿using MeshData;
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

//        private static string[] sampleFormulas =
//        {
//@"z < (sin(6*x)+sin(6*y))/4",
////------------------------------------------------------------
//@"(
//	x^2 + y^2 + z^2 < 0.75^2
//) ||
//(
//	x^2 + y^2 + (z-0.75)^2 < 0.25^2
//)",
////------------------------------------------------------------
//@"(
//    x^2 + y^2 + z^2 < 1.25 ^ 2
//) &&
//x^2+y^2 > 0.5^2 &&
//x^2+z^2 > 0.5^2 &&
//z^2+y^2 > 0.5^2",
////------------------------------------------------------------
//@"(
//    x^2 + y^2 + z^2 < 1 &&
//    x + y < z^2
//) ||
//(
//    x^2 + z^2 < 0.25
//)",
////------------------------------------------------------------
//@"(
//	x^2+y^2+z^2<1 &&
//	z<0+(cos(5*x)+cos(5*y))/3
//) ||
//x^2+y^2+(z-.8)^2<0.25",
////------------------------------------------------------------
//@"abs(x)<0.25 ||
//0.5*x^2+y^2+z^2<0.5",
////------------------------------------------------------------
//@"8*x^2+8*y^2+z^2<1 ||
//8*x^2+y^2+8*z^2<1 ||
//x^2+8*y^2+8*z^2<1",
////------------------------------------------------------------
//@"16*x^2+16*y^2+z^2<1 ||
//90*x^6+y^2+90*z^6<1 ||
//x^2+8*y^2+8*z^2<1",
////------------------------------------------------------------
//@"[advanced]
//public IneqMesh GetIneqMesh()
//{
//    return new IneqMesh
//    {
//        X0 = -1.5,
//        Y0 = -1.5,
//        Z0 = -1,
//        X1 = 1.5,
//        Y1 = 1.5,
//        Z1 = 5,
//        D = 0.15d,
//        Boxed = true,
//        IneqTree = IneqTreeParser.FromFormula(@""
//        (x-0.25*cos(4*z))^2+(y-0.25*sin(4*z))^2<0.85 &&
//        (x-0.25*cos(4*z))^2+(y-0.25*sin(4*z))^2>0.3
//        "")
//    };
//}
//",
////------------------------------------------------------------
//@"[advanced]
//public IneqMesh GetIneqMesh()
//{
//    return new IneqMesh
//    {
//        X0 = -1.0,
//        Y0 = -1.2,
//        Z0 = -1.3,
//        X1 = 1.0,
//        Y1 = 1.2,
//        Z1 = 1.3,
//        D = 0.1d,
//        Boxed = false,
//        IneqTree = IneqTreeParser.FromFormula(@""
//        (
//            (2*x^2+y^2+z^2-1)^3<(0.1*x^2+y^2)*z^3 &&
//            x<0.25
//        ) ||
//        (x^2+y^2+(z-0.2)^2<0.25)
//        "")
//    };
//}
//",
////------------------------------------------------------------
//@"[advanced]
//private IneqTree Balls(int count, double R, double r, double p)
//{
//    var res = new IneqTree((x, y, z) => 1);

//    for (int i = 0; i < count; i++)
//    {
//        double x0 = R * Math.Cos(i * 2 * Math.PI / count);
//        double y0 = R * Math.Sin(i * 2 * Math.PI / count);

//        res = res | ((x, y, z) =>  
//                Math.Pow(Math.Abs(x - x0),p) + 
//				Math.Pow(Math.Abs(y - y0),p) +
//				Math.Pow(Math.Abs(z),p) - 
//                Math.Pow(r,p));
//    }

//    return res;
//}

//public IneqMesh GetIneqMesh()
//{
//	var p = 5d;
//  	var ball = new IneqTree((x, y, z) => Math.Pow(Math.Abs(x),p) + Math.Pow(Math.Abs(y),p)  + Math.Pow(Math.Abs(z),p) - Math.Pow(0.7d,p));
  
//    return new IneqMesh
//    {
//        X0 = -1.0,
//        Y0 = -1.0,
//        Z0 = -1.0,
//        X1 = 1.0,
//        Y1 = 1.0,
//        Z1 = 1.0,
//        D = 0.1d,
//        Boxed = false,
//        IneqTree =
//                (
//                    ball &
//                   ((x, y, z) => -x*x - y*y + 0.05d)
//                ) |
//                Balls(8, 0.8d, 0.2d, p)
//    };
//}
//"
//        };

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

        //public int SampleIndex { get; set; }
        public int SampleUFuncIndex { get; set; }

        public static IneqMeshViewModel DefaultModel(int sample = -1, int sampleU = -1)
        {
            /*if (sample == -1)
            {
                sample = new Random().Next(sampleFormulas.Length);
            }
            else
            {
                sample = sample % sampleFormulas.Length;
            }*/

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
                Formula = "x^2 + y^2 + z^2 < 1",//sampleFormulas[sample],
                //SampleIndex = sample,
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
            if (Formula.StartsWith("[advancedX]"))
            {
                IneqMesh = DynamicCodeExecutor.Execute(Formula.Replace("[advancedX]", ""));
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

            if (!String.IsNullOrEmpty(Name))
            {
                string path = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/Samples"), Name + ".json");
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
                    Formula
                }, Formatting.None, new JsonSerializerSettings
                {
                    ContractResolver = new CamelCasePropertyNamesContractResolver(),
                    DateFormatHandling = DateFormatHandling.IsoDateFormat,
                    NullValueHandling = NullValueHandling.Ignore
                }));
            }
            //System.IO.File.WriteAllText(path, ineqMeshViewModel.PLY);
        }
    }
}