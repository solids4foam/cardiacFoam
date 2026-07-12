import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
MODEL = ROOT / "src/activeTensionModels/LandNiedererBatched/LandNiedererBatched.C"
HEADER = ROOT / "src/activeTensionModels/LandNiederer/LandNiederer_2017.H"


def test_batched_model_overrides_and_uses_hot_path_for_preconditioning():
    source = MODEL.read_text()
    declaration = MODEL.with_suffix(".H").read_text()

    assert "preconditionToRestingState" in declaration
    body = source.split("LandNiedererBatched::preconditionToRestingState", 1)[1]
    body = body.split("LandNiedererBatched::calculateTension", 1)[0]
    assert 'lookupOrDefault<scalar>("preconditioningTime", 1000.0)' in body
    assert "backend.evaluateScratchAtTime" in body
    assert "backend.applyExplicitEulerStep" in body
    assert "core_.scatterCellState" in body
    assert "prevLambda_ = 1.0" in body
    assert "lambdaRate_ = 0.0" in body


def test_generated_rhs_converges_with_default_batched_preconditioning_step(tmp_path):
    compiler = shutil.which("c++")
    if compiler is None:
        pytest.skip("C++ compiler unavailable")

    program = tmp_path / "precondition.cpp"
    executable = tmp_path / "precondition"
    program.write_text(
        f'''#include <algorithm>
#include <cmath>
#include <iostream>
#include "{HEADER}"

int main()
{{
    double constants[NUM_CONSTANTS];
    double rates[NUM_STATES];
    double ratesFine[NUM_STATES];
    double states[NUM_STATES];
    double statesFine[NUM_STATES];
    double algebraic[NUM_ALGEBRAIC] = {{}};
    LandNiederer2017initConsts(constants, rates, states);
    LandNiederer2017initConsts(constants, ratesFine, statesFine);

    constexpr double restingCai = 0.0002;
    constexpr double step = 0.01;
    constexpr int nSteps = 100000;
    for (int stepI = 0; stepI < nSteps; ++stepI)
    {{
        algebraic[AV_Cai] = restingCai;
        algebraic[AV_lambda] = 1.0;
        algebraic[AV_lambda_rate] = 0.0;
        LandNiederer2017computeVariables
        (
            stepI*step, constants, rates, states, algebraic
        );
        for (int stateI = 0; stateI < NUM_STATES; ++stateI)
        {{
            states[stateI] += step*rates[stateI];
        }}
        for (int halfI = 0; halfI < 2; ++halfI)
        {{
            LandNiederer2017computeVariables
            (
                stepI*step + halfI*step*0.5,
                constants, ratesFine, statesFine, algebraic
            );
            for (int stateI = 0; stateI < NUM_STATES; ++stateI)
            {{
                statesFine[stateI] += step*0.5*ratesFine[stateI];
            }}
        }}
    }}

    algebraic[AV_Cai] = restingCai;
    algebraic[AV_lambda] = 1.0;
    algebraic[AV_lambda_rate] = 0.0;
    LandNiederer2017computeVariables
    (
        1000.0, constants, rates, states, algebraic
    );
    double maxRate = 0.0;
    for (double rate : rates) maxRate = std::max(maxRate, std::abs(rate));
    double maxStepDifference = 0.0;
    for (int stateI = 0; stateI < NUM_STATES; ++stateI)
    {{
        maxStepDifference = std::max
        (
            maxStepDifference,
            std::abs(states[stateI] - statesFine[stateI])
        );
    }}
    std::cout << maxRate << " " << maxStepDifference << " "
              << algebraic[AV_Ta] << "\\n";
}}
'''
    )
    subprocess.run(
        [compiler, "-std=c++17", "-O2", str(program), "-o", str(executable)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [str(executable)], check=True, capture_output=True, text=True
    )
    max_rate, max_step_difference, tension = map(float, result.stdout.split())
    assert max_rate < 1e-8
    assert max_step_difference < 1e-8
    assert tension > 0.0
