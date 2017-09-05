# Reliable-Re-detection-for-Long-term-Tracking
Code for the re-detection framework proposed in the paper **Reliable Re-detection for long-term tracking**

The re-detection tracker is based on the Staple, and the deep re-detection tracker takes HCF as baseline. 

For the re-detection tracker, just start Matlab and run the **runTracker.m**. To run the latter, deep re-detection tracker, please download the VGG-19 and compile the Matconvnet following the description in its **README**.  

This code includes a quite general framework for long-term tracking. You can easily incorporate your own DCF based trackers.

About speed and performance:

Since our framework re-utilizes the baseline tracker for re-detection, the speed is closely related to the baseline method as well as the difficulty degree of the video.
With only hand-crafted features, the re-detection tracker achieves more than 40 FPS on a single CPU.

Tracking performance may vary slightly (less than 1%) on different machines.
This is due to small numerical effects which can accumulate over time (actually, all trackers are affected by this) and the random errors caused by particle filter.
After increasing the particle number, the performance will improve slightly and be more stable. 

# Contact
If you find any questions, please feel free to contact wn6149@mail.ustc.edu.cn

More details will be introduced latter. 
