##Topic: Offset-free Nonlinear Model Predictive Control: A comparison of different formulation


## Statement of Research Problem
The research problem at the core of this study is the ongoing challenge of achieving offset-free Nonlinear Model Predictive Control (NMPC). Despite the proven effectiveness of NMPC in managing complex nonlinear systems, 
one of its persistent issues is the presence of steady-state offsets, which lead to suboptimal control performance and deviations from desired setpoints. 
This issue is particularly problematic when there are model-plant mismatches or unmeasured disturbances.
To address this challenge, this research focuses on comparing two prominent state estimation techniques—Extended Kalman Filter (EKF) and Moving Horizon Estimation (MHE)—which are both widely used to improve the accuracy of state estimation and enhance the robustness of NMPC. 

The objective of this study is to determine which formulation offers a more efficient and reliable solution for achieving offset-free control in nonlinear systems. The research will involve developing simulation models of nonlinear processes(Continuous Fermenter, Continuous Stirred Tank Reactor, and Quadruple Tank System), 
implementing NMPC controllers with EKF and MHE, and conducting a detailed comparative analysis of their performance. This comparison aims to identify the most effective and robust approach for offset-free NMPC, potentially contributing to advancements in the control of complex nonlinear systems.



## This repository contains MATLAB and Simulink code for simulating a **Offset-free Nonlinear Model Predictive Control (NMPC)** approach applied to a **Nonlinear processes**. The simulation response show how the system response to disturbance and process model mismatch.

## A README.m FILE IS ATTACHED TO EACH FOLDER TO HELP IN UNDERSTANDING AND SIMULATION OF THE CODE.

---

## 📁 Directory Structure

```plaintext
CON_FER/
README.m
├── MHE/
│   └── state/
│       ├── NMPC_conFer.m
│       ├── SIM_continuousFermenter.slx
│       ├── conFerStateFcnCT.m
│       ├── conFerOutputFcn.m
│       ├── MHE_compute.m
│       ├── plots.m
│       ├── response.jpg
│       ├── estimated_state.jpg
│       ├── state_dist.jpg
│       └── ...
```

> Ensure you check the README.m file in each folder before runing any simulation.

---


## 🛠 Requirements

- **MATLAB** (Recommended: R2023b or later)
- **Simulink**
- **Model Predictive Control Toolbox**

---

## 📄 License

This project is intended for educational and research purposes only.

---

## 📬 Contact

**Author:** Hammed Akeeb  
**Email:** hammedxyz@gmail.com
