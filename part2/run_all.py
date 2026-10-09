"""Run all four experiments; figures go to figures/, summaries to results/."""

import exp1_amplitude
import exp2_energy
import exp3_damped_accuracy
import exp4_damped_nonlinear

if __name__ == "__main__":
    for module in (exp1_amplitude, exp2_energy, exp3_damped_accuracy, exp4_damped_nonlinear):
        print(f"\n### running {module.__name__}\n")
        module.main()
