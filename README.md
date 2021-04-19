# Robust MPC for Linear Systems with Parametric and Additive Uncertainty: A Novel Constraint Tightening Approach 

These set of codes replicate the Region of Attraction (ROA) and the computation times presented in the draft https://arxiv.org/abs/2007.00930. The yellow set denotes the ROA of our MPC, which is compared to the corresponding ROAs from Tube MPC (grey set, www.sciencedirect.com/science/article/pii/S0005109803002838?via%3Dihub, Section 5), and System level Synthesis based Constrained LQR (blue set, https://arxiv.org/pdf/1809.10121.pdf, Section 2.3).  

![roa](https://user-images.githubusercontent.com/12418616/115236563-571b9800-a0d0-11eb-8509-ce53fdce4aec.png)
The ROA of our approach is 12x bigger than the ROA from the Constrained LQR. Although the ROA of the Tube MPC is comparable by volume, our algorithm (denoted by Algorithm 1), is computationally more efficient than the Tube MPC for all considered horizons. This is seen in the table below. 

![tab](https://user-images.githubusercontent.com/12418616/115236885-c1343d00-a0d0-11eb-90a2-c0c96e110da8.png)
