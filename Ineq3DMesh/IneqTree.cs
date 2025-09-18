using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MeshData
{
    using FuncXYZ = Func<double, double, double, double>;
    public delegate void Transformer(ref double x, ref double y, ref double z);

    public class IneqTree
    {
        private IneqNode root = null;
        private List<FuncXYZ> expressionList = null;

        public List<FuncXYZ> ExpressionList
        {
            get 
            {
                if (expressionList == null)
                {
                    expressionList = new List<FuncXYZ>();
                    AddToexpressionList(root);
                }

                return expressionList; 
            }
        }

        public IneqNode Root
        {
            get { return root; }
        }

        public IneqTree()
        { 
        }

        public IneqTree(NodeType nodeType, IneqTree left, IneqTree right)
        {
            if (left.root == null)
            {
                root = right.root;
            }
            else if (right.root == null)
            {
                root = left.root;
            }
            else
            {
                root = new IneqNode();
                root.Expression = null;
                root.NodeType = nodeType;
                root.Left = left.root;
                root.Right = right.root;
            }
        }

        public IneqTree(FuncXYZ expression)
        {
            root = new IneqNode();
            root.Expression = expression;
            root.NodeType = NodeType.NodeExpression;
            root.Left = null;
            root.Right = null;
        }

        private void AddToexpressionList(IneqNode node)
        {
            if (node == null)
                return;

            if (node.NodeType == NodeType.NodeExpression)
            {

                /*int i = 0;
                while (i < expressionList.Count && !Equal(node.Expression, expressionList[i]))
                {
                    i++;
                }

                if (i == expressionList.Count)*/
                {
                    expressionList.Add(node.Expression);
                    //node.Expression = null;
                    node.ExpressionIndex = expressionList.Count - 1;
                }
                /*else
                {
                    node.ExpressionIndex = i;
                    node.Expression = expressionList[i];
                }*/
            }

            AddToexpressionList(node.Left);
            AddToexpressionList(node.Right);
        }

        /*private bool Equal(FuncXYZ f1, FuncXYZ f2)
        {
            Random rand = new Random();
            int p = 0;
            while (p < 100)
            {
                double x = rand.NextDouble() * 2 - 1;
                double y = rand.NextDouble() * 2 - 1;
                double z = rand.NextDouble() * 2 - 1;

                if (Math.Abs(f1(x, y, z) - f2(x, y, z)) > 1e-10)
                {
                    return false;
                }

                p++;
            }
            return true;
        }*/

        public IneqTree Not()
        {
            if (root != null)
            {
                root.Not();
            }

            return this;
        }

        public IneqTree Transform(Transformer t)
        {
            if (root != null)
            {
                root.Transform(t);
            }

            return this;
        }

        public IneqTree ToFunc(double eps)
        {
            if (root != null)
            {
                return new IneqTree(root.ToFunc(eps));
            }
            else
            { 
                return null;
            }
        }

        public static IneqTree operator !(IneqTree ineqTree)
        {
            return ineqTree.Not();
        }

        public static IneqTree operator &(IneqTree left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator |(IneqTree left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static IneqTree operator &(IneqTree left, FuncXYZ right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator &(FuncXYZ left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeAnd, left, right);
        }

        public static IneqTree operator |(IneqTree left, FuncXYZ right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static IneqTree operator |(FuncXYZ left, IneqTree right)
        {
            return new IneqTree(NodeType.NodeOr, left, right);
        }

        public static implicit operator IneqTree(FuncXYZ expression)
        {
            return new IneqTree(expression);
        }
        
        public enum NodeType
        {
            NodeExpression,
            NodeOr,
            NodeAnd
        }

        public class IneqNode
        {
            private NodeType nodeType;
            private FuncXYZ expression = null;
            private int expressionIndex = -1;
            private IneqNode left = null;
            private IneqNode right = null;

            public NodeType NodeType
            {
                get { return nodeType; }
                set { nodeType = value; }
            }

            public FuncXYZ Expression
            {
                get { return expression; }
                set { expression = value; }
            }

            public int ExpressionIndex
            {
                get { return expressionIndex; }
                set { expressionIndex = value; }
            }

            public IneqNode Left
            {
                get { return left; }
                set { left = value; }
            }

            public IneqNode Right
            {
                get { return right; }
                set { right = value; }
            }

            public IEnumerable<int> ExpressionIndexes
            {
                get
                {
                    if (nodeType == IneqTree.NodeType.NodeExpression)
                        return new int[] { expressionIndex };
                    else
                        return left.ExpressionIndexes.Union(right.ExpressionIndexes);
                }
            }

            public int SubNodeCount
            {
                get
                {
                    if (nodeType == IneqTree.NodeType.NodeExpression)
                        return 0;
                    else
                        return left.SubNodeCount + right.SubNodeCount;
                
                }
            }

            public int DomainCount
            {
                get
                {
                    if (nodeType == IneqTree.NodeType.NodeExpression)
                        return 1;
                    else
                        return left.DomainCount + right.DomainCount + 1;

                }
            }

            internal void Not()
            {
                if (NodeType == NodeType.NodeExpression)
                {
                    var oldEx = Expression;
                    Expression = ((x, y, z) => - oldEx(x, y, z));
                }
                else if (NodeType == NodeType.NodeAnd)
                {
                    NodeType = NodeType.NodeOr;
                    left.Not();  right.Not();
                }
                else// NodeType == NodeType.NodeOr
                {
                    NodeType = NodeType.NodeAnd;
                    left.Not(); right.Not();                    
                }
            }

            internal void Transform(Transformer t)
            {
                if (NodeType == NodeType.NodeExpression)
                {
                    var oldEx = Expression;
                    Expression = ((x, y, z) => 
                    {
                        t(ref x, ref y, ref z);
                        return oldEx(x, y, z);
                    });
                }
                else
                {
                    left.Transform(t); right.Transform(t);
                }
            }

            internal FuncXYZ ToFunc(double eps)
            {
                if (NodeType == NodeType.NodeExpression)
                {
                    return Expression;
                }

                var l = left.ToFunc(eps);
                var r = right.ToFunc(eps);

                if (NodeType == NodeType.NodeAnd)
                {
                    if (eps > 0)
                    {
                        return ((x, y, z) => 
                        {
                            var a = l(x, y, z);
                            var b = r(x, y, z);
                            return (a + b + SmoothAbs(a - b,eps)) / 2.0d;
                        });
                    }
                    else
                    {
                        return ((x, y, z) => { return Math.Max(l(x, y, z), r(x, y, z)); });
                    }
                }
                else 
                {
                    if (eps > 0)
                    {
                        return ((x, y, z) =>
                        {
                            var a = l(x, y, z);
                            var b = r(x, y, z);
                            return (a + b - SmoothAbs(a - b, eps )) / 2.0d;
                        });

                    }
                    else
                    {
                        return ((x, y, z) => { return Math.Min(l(x, y, z), r(x, y, z)); });
                    }
                }
            }

            static double Bump0 = Bump(0.0);
            public static double SmoothAbs(double x, double epsilon0 = 0.03d)
            {
                
                double R = 10 * epsilon0;

                double ax = Math.Abs(x);
                if (ax >= R)
                {
                    return ax; // exact outside transition zone
                }

                // Normalized t in [0,1]
                double t = ax / R;

                // Cubic Hermite polynomial: ψ(t) = 1 - t - t^2 + t^3
                //double psi = 1 - t * t + t * t * t;

                // ε(|x|), smoothly decaying from ε0 at 0 to 0 at R
                //double eps = epsilon0 * psi;

                // Normalized bump: b(t)/b(0)
                double bump = Bump(t) / Bump0;

                double eps = epsilon0 * bump;

                // Smooth approximation
                return Math.Sqrt(x * x + eps * eps);
            }

            private static double Bump(double t)
            {
                if (Math.Abs(t) >= 1.0) return 0.0;
                double denom = 1.0 - t * t;
                return Math.Exp(-1.0 / denom);
            }
        }
    }
}
