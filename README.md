# Robust MPC for Linear Systems with Parametric and Additive Uncertainty: A Novel Constraint Tightening Approach 

These set of codes replicate the Region of Attraction (ROA) and the computation times presented in the draft https://arxiv.org/abs/2007.00930. The yellow set denotes the ROA of our MPC, which is compared to the corresponding ROAs from Tube MPC (grey set, www.sciencedirect.com/science/article/pii/S0005109803002838?via%3Dihub, Section 5), and System level Synthesis based Constrained LQR (blue set, https://arxiv.org/pdf/1809.10121.pdf, Section 2.3).  

![roa](https://user-images.githubusercontent.com/12418616/116338475-a593fb00-a790-11eb-9f2a-a8abcc0d128e.png)

The ROA of our approach (denoted by Algorithm 1) is ~12x bigger than the ROA from the Constrained LQR. Although the ROA of the Tube MPC is only about 4% smaller in volume, our algorithm, is computationally more efficient than the Tube MPC for all considered horizons. This is seen in the table below. The offline times required in our algorithm are for pre-computation of a set of bounds before control design. Notice that our offline + online combined times are still faster than the online times needed for the Tube MPC.

![tim](https://user-images.githubusercontent.com/12418616/116338519-b3498080-a790-11eb-9e7d-9aa28852b418.png)
