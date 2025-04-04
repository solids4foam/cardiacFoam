import numpy as np

from cardiac_benchmark import activation_model
from cardiac_benchmark.problem import LVProblem


def test_active_model_default_parameters():
    problem_parameters = LVProblem.default_parameters()
    active_parameters = activation_model.default_parameters()
    dt = float(problem_parameters["dt"])
    time = np.arange(dt, 1, dt)
    t_eval = time - float(problem_parameters["alpha_f"]) * dt
    tau = activation_model.activation_function(
        t_span=(0, 1),
        t_eval=t_eval,
        parameters=active_parameters,
    )
    # Different value than in the paper, but it depends on dt, alpha_f and integration method
    assert np.isclose(tau.max(), 118068.74829471917)

    # import matplotlib.pyplot as plt

    # plt.plot(time, tau)
    # plt.savefig("tau.png")
