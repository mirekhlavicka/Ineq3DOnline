using Jace;
using MeshData;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading;
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
            error = false;
            IneqTree result = null;

            try
            {
                result = BuildIneqTree(formulaText
                    .Replace("&&", "&")
                    .Replace("||", "|")
                    .Replace(" ", "")
                    .Replace("\n", "")
                    .Replace("\r", "")
                    .Replace("\t", "")
                    .Replace(",", "."));
            }
            catch (Exception exc)
            {
                if (!error)
                {
                    throw exc;
                }
            }

            if (error)
            {
                throw new Exception("Syntax error in expression");
            }

            return result;
        }

        private static bool error = false;
        private static IneqTree BuildIneqTree(string vyraz)
        {
            int i, pocet_zavorek = 0, pozice_operatoru = -1;
            bool nasel_se_and = false;
            bool nasel_se_or = false;
            bool nasla_se_ner = false;
            char z;  // zkoumany znak retezce
            IneqTree vysledek = null, node_levy = null, node_pravy = null;

            error = error || (vyraz.Length == 0);

            if (!error)
            {
                i = 0;
                while ((i < vyraz.Length) && (!nasel_se_or))
                {
                    z = vyraz[i];
                    if (z == '(') pocet_zavorek++;
                    else if (z == ')') pocet_zavorek--;
                    else if (pocet_zavorek == 0)
                    {
                        //budu zkoumat, jestli se nenasel zajimavy logicky operator nebo nerovnost
                        if (z == '|')
                        {
                            // nasel se OR, zapamatuji si kde a cyklus skonci
                            nasel_se_or = true;
                            pozice_operatoru = i;
                        }
                        else if (!nasel_se_and)   // jestli se uz nejake AND naslo, tak uz me zajima jenom pripadne nasledne OR
                        {
                            if (z == '&')
                            {
                                // nasel se prvni and, zapamatuji si kde a budu hledat dal, jestli tam treba neni or
                                pozice_operatoru = i;
                                nasel_se_and = true;
                            }
                            else if (!nasla_se_ner) // naslo-li se nejake < nebo >, hledam uz jenom logicke operatory
                            {
                                if ((z == '<') || (z == '>'))
                                {
                                    pozice_operatoru = i;
                                    nasla_se_ner = true;
                                }
                            }  // if (!nasla_se_ner)
                        }  // if (!nasel_se_and)
                    }  // if (pocet_zavorek==0)
                    i++;
                } // while

                if (pocet_zavorek != 0)  // je to blbe, nesouhlasi pocet zavorek
                    error = true;
                else if (nasel_se_or)
                {
                    // vytvorim uzel a rozdelim retezec a pustim na jeho casti stromovac
                    string levy = vyraz.Substring(0, pozice_operatoru);
                    string pravy = vyraz.Substring(pozice_operatoru + 1, vyraz.Length - pozice_operatoru - 1);
                    node_levy = BuildIneqTree(levy);
                    node_pravy = BuildIneqTree(pravy);
                    vysledek = new IneqTree(IneqTree.NodeType.NodeOr, node_levy, node_pravy);
                }
                else if (nasel_se_and)
                {
                    // dtto
                    string levy = vyraz.Substring(0, pozice_operatoru);
                    string pravy = vyraz.Substring(pozice_operatoru + 1, vyraz.Length - pozice_operatoru - 1);
                    node_levy = BuildIneqTree(levy); node_pravy = BuildIneqTree(pravy);
                    vysledek = new IneqTree(IneqTree.NodeType.NodeAnd, node_levy, node_pravy);
                }
                else if (nasla_se_ner)
                {
                    //  upravim nerovnost, ukoncim vetev, ulozim retezec
                    PrepareFormula(ref vyraz, pozice_operatoru);

                    Func<double, double, double, double> formula = (Func<double, double, double, double>)engine.Formula(vyraz)
                        .Parameter("x", DataType.FloatingPoint)
                        .Parameter("y", DataType.FloatingPoint)
                        .Parameter("z", DataType.FloatingPoint)
                        .Result(DataType.FloatingPoint)
                        .Build();

                    vysledek = new IneqTree(formula);
                }
                else // vsechny zavorky jsou zavrene, ale nenasel se operator
                {
                    // mohou nastat dve moznosti: jsou tam zavorky navic nebo je to blbe
                    if ((vyraz[0] == '(') && (vyraz[vyraz.Length - 1] == ')'))
                    {
                        string kratsi = vyraz.Substring(1, vyraz.Length - 2);
                        vysledek = BuildIneqTree(kratsi);
                    }
                    else error = true;
                } // nenasel se operator

            } // if ( !error )

            return vysledek;

        }  // funkce

        private static void PrepareFormula(ref string vyraz, int pozice_nerovnosti)
        {
            string leva = vyraz.Substring(0, pozice_nerovnosti);
            string prava = vyraz.Substring(pozice_nerovnosti + 1, vyraz.Length - pozice_nerovnosti - 1);
            if ( vyraz[pozice_nerovnosti] == '<' )
            {
                if (prava == "0") vyraz = leva;   // je to neco oper 0
                else if (leva == "0") vyraz = prava;   // je to 0 oper neco
                else vyraz = leva + "-(" + prava + ")";
            }
            else // je to neco > neco
            {
                if (leva == "0") vyraz = prava;   //  je to 0 > neco
                else vyraz = prava + "-(" + leva + ")";
            }
        }
    }
}