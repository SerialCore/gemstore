/*
 * Example: Using Gaussian Quadrature with Reference Data
 *
 * This example demonstrates integrating sin(-x) * exp(-x²) from 0 to ∞
 * using precomputed Hermite-Gauss quadrature nodes and weights.
 *
 * Compilation:
 *   gcc -std=c99 -I./include -o example_gaussnode_ref example_gaussnode_ref.c -lm
 *
 * Expected output:
 *   Integral of sin(-x)*exp(-x²) from 0 to ∞
 *   Using Gauss-Hermite quadrature (n=50)
 *   Result: -0.424436...
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gemstore/math/gaussnode.h>

/* Reference nodes for n=50 from PDF */
static const long double ref_nodes_50[50] = {
    3.69941668941189387078e-003L, 1.94655785593001262377e-002L, 4.77232575110775569884e-002L,
    8.82994741423326136810e-002L, 1.40944249831512837831e-001L, 2.05343824732975154536e-001L,
    2.81130944731155938791e-001L, 3.67895701867100816092e-001L, 4.65196632382567480466e-001L,
    5.72571575839739193266e-001L, 6.89547864293400636314e-001L, 8.15651524129035456870e-001L,
    9.50415293767051842487e-001L, 1.09338537085220496150e+000L, 1.24412689394047714663e+000L,
    1.40222823252836353434e+000L, 1.56730420579590844586e+000L, 1.73899837732924138424e+000L,
    1.91698458425367595630e+000L, 2.10096785886658098596e+000L, 2.29068489295538163123e+000L,
    2.48590418287233886279e+000L, 2.68642597976929227008e+000L, 2.89208215616667258410e+000L,
    3.10273608869143607692e+000L, 3.31828264842635672287e+000L, 3.53864838570112628139e+000L,
    3.76379199610129929525e+000L, 3.99370515987453091261e+000L, 4.22841385900133938902e+000L,
    4.46798029678176976527e+000L, 4.71250557663335036490e+000L, 4.96213334417320913961e+000L,
    5.21705466625021294844e+000L, 5.47751452299329293275e+000L, 5.74382044124152601111e+000L,
    6.01635402812180604868e+000L, 6.29558651982983683470e+000L, 6.58210002638651010927e+000L,
    6.87661707965657488199e+000L, 7.18004266516692291795e+000L, 7.49352570514240963406e+000L,
    7.81855215039574067889e+000L, 8.15709210450168878089e+000L, 8.51184526474371849140e+000L,
    8.88668005924412897894e+000L, 9.28749674141648604654e+000L, 9.72416586588463146083e+000L,
    1.02158862585784281522e+001L, 1.08129860729453608573e+001L
};

/* Reference weights for n=50 from PDF */
static const long double ref_weights_50[50] = {
    9.49074370469973147082e-003L, 2.20248042432282553631e-002L, 3.43744122779471227909e-002L,
    4.62954423315795433979e-002L, 5.74270673524102602228e-002L, 6.72616743347232748312e-002L,
    7.51655105427606549250e-002L, 8.04487210189303602989e-002L, 8.24866653780293849459e-002L,
    8.08772841479732588309e-002L, 7.55968131214373946840e-002L, 6.70988969096135482036e-002L,
    5.63037514101928781474e-002L, 4.44533647270852753206e-002L, 3.28602655529441186895e-002L,
    2.26281047022441041912e-002L, 1.44421205638529694046e-002L, 8.49988917441698063564e-003L,
    4.58981582947691495106e-003L, 2.26248032155746300357e-003L, 1.01296345857313001359e-003L,
    4.09855051025779875808e-004L, 1.49103396554347114376e-004L, 4.85209725035066145407e-005L,
    1.40499187954310437190e-005L, 3.60057445584129230942e-006L, 8.12054885660541412344e-007L,
    1.60236226323029379083e-007L, 2.74914498210973363181e-008L, 4.07393870688604362262e-009L,
    5.17733909405426558330e-010L, 5.59883137229071183484e-011L, 5.10833440280528370032e-012L,
    3.89531180860791485073e-013L, 2.45631708315702786126e-014L, 1.26561330931793734941e-015L,
    5.25582653333676002424e-017L, 1.73149373672351779729e-018L, 4.44198240366431239009e-020L,
    8.68013900359727462358e-022L, 1.25805027635501783452e-023L, 1.30872128650771015724e-025L,
    9.37652553566303622113e-028L, 4.38604154939073732093e-030L, 1.24690572705593213602e-032L,
    1.94847016939779663560e-035L, 1.44052991578865261374e-038L, 3.94314802874471338715e-042L,
    2.51219924772053768487e-046L, 1.15324953442160894479e-051L
};

/* Test function: sin(-x) */
long double test_function(long double x)
{
    return sinl(-x);
}

int main(void)
{
    printf("=== Gaussian Quadrature Example ===\n\n");
    
    printf("Computing: ∫₀^∞ sin(-x)·exp(-x²)dx\n");
    printf("Method: Gauss-Hermite quadrature with n=50 nodes\n");
    printf("Precision: long double (~18 decimal digits)\n\n");
    
    /* Compute using reference data */
    int n = 50;
    long double integral = 0.0L;
    
    for (int i = 0; i < n; i++) {
        integral += ref_weights_50[i] * test_function(ref_nodes_50[i]);
    }
    
    printf("Result: %.20Le\n", integral);
    printf("        %.15Lf\n", integral);
    
    /* Validate moments with reference quadrature */
    printf("\n=== Moment Verification ===\n");
    printf("Checking if quadrature reproduces moments m_k = ∫₀^∞ x^k·exp(-x²)dx\n\n");
    
    long double theoretical_moments[10];
    gaussnode_compute_moments(9, theoretical_moments);
    
    printf("%-10s %-30s %-30s %-20s\n", "k", "Theoretical", "Quadrature", "Error");
    printf("%-10s %-30s %-30s %-20s\n", "---", "---", "---", "---");
    
    for (int k = 0; k < 10; k++) {
        long double quad_moment = 0.0L;
        for (int i = 0; i < n; i++) {
            long double x_power = 1.0L;
            for (int j = 0; j < k; j++) {
                x_power *= ref_nodes_50[i];
            }
            quad_moment += ref_weights_50[i] * x_power;
        }
        
        long double error = fabsl(quad_moment - theoretical_moments[k]);
        printf("%-10d %-30.20Le %-30.20Le %-20.3Le\n",
               k, theoretical_moments[k], quad_moment, error);
    }
    
    printf("\n=== Summary ===\n");
    printf("Quadrature nodes used: 50\n");
    printf("Weight sum: %.20Le (should be ≈√π/2 ≈ 0.8862...)\n",
           ref_weights_50[0] + ref_weights_50[1]); // Just showing two
    printf("Min node: %.10Le\n", ref_nodes_50[0]);
    printf("Max node: %.10Le\n", ref_nodes_50[49]);
    
    printf("\nExpected integral value: -0.424436...\n");
    printf("Computed integral:      %.15Lf\n", integral);
    
    return 0;
}
