using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace VisualPSV
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        //Get normmal distibution using a Box-Muller transform
        private double RandomNormal(double mean, double std) 
        {
            Random rand = new Random(); //reuse this if you are generating many
            double u1 = 1.0 - rand.NextDouble(); //uniform(0,1] random doubles
            double u2 = 1.0 - rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            double randNormal =
                         mean + std * randStdNormal; //random normal(mean,std^2)

            return randNormal;
        }

        //Returns sigma, corresponding to hamonic l
        private double SL(int l, double alpha = 2.7)
        {

            double a2 = Math.Pow(alpha, 2);

            double ca = 0.547;

            double s_l = Math.Sqrt(a2*Math.Pow(ca,2*l)/((l+1)*(2*l+1)));
            return s_l;
        }

        // This functions returns h and g coeffs
        private List<double> getCoeffs(int terms, double G1, double G2, double G3, double beta, double afact)
        {
            List<double> gh = new List<double>();

            double alpha = G1/afact;

            double s1 = SL(1, alpha);
            
            double s10 = s1 * beta;

            gh.Add(RandomNormal(G1,s1));
            gh.Add(RandomNormal(0,s1));
            gh.Add(RandomNormal(0, s1));

            double s = 0;
            double o = 0;

            for (int l = 2; l < (terms+1);l++)
            {
                for (int m = 0; m < l+1; m++)
                {
                    o = 0;

                    if(l == 2 && m == 0)
                    {o = G2;}

                    if (l == 3 && m == 0)
                    { o = G3; }

                    s = SL(l,alpha);

                    if ((l - m) % 2 == 1)
                    { s = s * beta; }

                    gh.Add(RandomNormal(o,s));

                    if(m != 0)
                    {gh.Add(RandomNormal(0, s));}
                }
            }

            return gh;
        }

    }
}
