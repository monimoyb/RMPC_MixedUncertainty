# Robust MPC for Linear Systems with Parametric and Additive Uncertainty: A Novel Constraint Tightening Approach 

These set of codes replicate the Region of Attraction (ROA) and the computation times presented in the draft https://arxiv.org/abs/2007.00930. The yellow set denotes the ROA of our MPC, which is compared to the corresponding ROAs from Tube MPC (grey set, www.sciencedirect.com/science/article/pii/S0005109803002838?via%3Dihub, Section 5), and System level Synthesis based Constrained LQR (blue set, https://arxiv.org/pdf/1809.10121.pdf, Section 2.3).  

![roa](https://user-images.githubusercontent.com/12418616/115521712-468a2f80-a240-11eb-9e5f-a9011b7233c3.png)

The ROA of our approach (denoted by Algorithm 1) is ~12x bigger than the ROA from the Constrained LQR. Although the ROA of the Tube MPC is only about 4% smaller in volume, our algorithm, is computationally more efficient than the Tube MPC for all considered horizons. This is seen in the table below. 

![time](https://user-images.githubusercontent.com/12418616/115521737-4e49d400-a240-11eb-9193-354944b86bc7.png)
