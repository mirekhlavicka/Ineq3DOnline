using Jace;
using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline
{
    public class IneqTreeParser //: IneqTree
    {
        //https://github.com/pieterderycke/Jace/wiki
        //http://www.codeproject.com/Articles/682589/Jace-NET-Just-another-calculation-engine-for-NET
        private static CalculationEngine engine = new CalculationEngine();

        public static IneqTree FromFormula(string formulaText)
        {
            if (formulaText.Contains("||"))
            {
                string[] tmp = formulaText.Split(new string[] { "||" }, StringSplitOptions.None);

                return new IneqTree(IneqTree.NodeType.NodeOr, FromFormula(tmp[0]), FromFormula(tmp[1]));
            }
            else if (formulaText.Contains("&&"))
            {
                string[] tmp = formulaText.Split(new string[] { "&&" }, StringSplitOptions.None);

                return new IneqTree(IneqTree.NodeType.NodeAnd, FromFormula(tmp[0]), FromFormula(tmp[1]));
            }
            else
            {
                Func<double, double, double, double> formula = (Func<double, double, double, double>)engine.Formula(formulaText.Replace("<", "-"))
                    .Parameter("x", DataType.FloatingPoint)
                    .Parameter("y", DataType.FloatingPoint)
                    .Parameter("z", DataType.FloatingPoint)
                    .Result(DataType.FloatingPoint)
                    .Build();

                return new IneqTree(formula);
            }
        }
    }
}