Here are the Time merching solution with KF Estimated one using model reduction for plate type of problem



<img width="1190" height="661" alt="Screenshot 2025-10-12 184337" src="https://github.com/user-attachments/assets/9a690343-bfec-4ee3-80e0-9a98ce6d0d8e" />

<img width="1104" height="622" alt="Screenshot 2025-10-12 184602" src="https://github.com/user-attachments/assets/cc03dc41-fb3f-4fa3-a7ff-3e5cfe790d77" />
<img width="1246" height="621" alt="Screenshot 2025-10-14 120518" src="https://github.com/user-attachments/assets/04463a57-4b4e-4ca8-bbd2-ff718e4214fc" />


After Impleneting the time merching solution for plate like structure in FEM ( using Neumark Beta Integration Scheme)  I get this results :

( In this plot I have used different different cases for weightage - i. biased towards measured value (0.1, 0.9),  ii. Balaced Approach(0.6,0.4)  iii. Focusing Fem Predicted Solution (0.9,0.1)
<img width="1104" height="622" alt="Screenshot 2025-10-12 184602" src="https://github.com/user-attachments/assets/e9c5a24b-b6f3-46c9-87e4-80cc9b279473" />

 While Plotting the Measured points and comparing the results with FEM predicted and KF estimated Solution
 (I know there are some mistakes while comparing the estimating the kalman filter feedback loop it will be easily debugged)

<img width="886" height="644" alt="Screenshot 2025-10-12 184313" src="https://github.com/user-attachments/assets/e8732971-9411-4acd-80eb-49d22912fad9" />

<img width="1035" height="670" alt="Screenshot 2025-10-12 184322" src="https://github.com/user-attachments/assets/548d94f5-4928-4139-88c1-db5af38724d7" />

<img width="1190" height="661" alt="Screenshot 2025-10-12 184337" src="https://github.com/user-attachments/assets/f2d97de9-7024-44f6-b00f-a1bd969cf9ae" />


<h/>


Results of  KF_PLATE(beamExtended).m
<img width="1175" height="649" alt="Screenshot 2025-10-08 170916" src="https://github.com/user-attachments/assets/cc2f98ae-77de-4a72-b32e-ff9c0f56aaa2" />
<img width="1146" height="677" alt="Screenshot 2025-10-07 203857" src="https://github.com/user-attachments/assets/9a5e2dc8-5912-46de-99ea-320ebbfb5a55" />
<img width="1031" height="647" alt="Screenshot 2025-10-08 182105" src="https://github.com/user-attachments/assets/e2cc374f-34f8-47ab-a8b3-40ac299c158b" />





While Starting With KF I have practiced do some handson the 1D FEM with its Simplest form ( Spring Mass System) Then I move on the other
Results form KF_Spring_noise_conrol.m

<img width="1162" height="562" alt="Screenshot 2025-10-03 094226" src="https://github.com/user-attachments/assets/63399c5f-dce9-4953-9624-1ed650535736" />
<img width="1209" height="578" alt="Screenshot 2025-10-03 093448" src="https://github.com/user-attachments/assets/1d84a36e-1d46-4dcc-ad16-306c01330588" />
<img width="834" height="366" alt="Screenshot 2025-10-03 094553" src="https://github.com/user-attachments/assets/1e146969-76ab-4ce7-bc70-e628a80e3231" />
<img width="1188" height="583" alt="Screenshot 2025-10-03 093834" src="https://github.com/user-attachments/assets/ca60bb94-31ee-4cbe-8d59-d78c5ef48618" />
