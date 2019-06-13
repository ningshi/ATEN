# ATEN
And/Or Tree Ensemble for inferring accurate Boolean network topology and dynamics 

Note: In the previous version, please note that set.seed() is not suitable for the Package parallel. Instead, we introduced function clusterSetRNGStream(cl,iseed) to make the results reproducible, more deatails see Pacakge parallel.

```
# Install devtools from CRAN
install.packages("devtools")

# Download the development version of ATEN from GitHub:
devtools::install_github("ningshi/ATEN")

# Load it to library
library("ATEN")
```


### Manual/Usage
--- 

<b>And/Or tree</b>: <br/>
We use list to represent an And/Or tree (i.e., a Boolean function). For instance, assuming we have a network consisting of 5 nodes in which the target node is x_1, and the Boolean function f we will infer is: f(x_1) = x1 || x2&&!x3 || !x2&&x5. Here, please note that we use time-series data, namely we have x_1(t+1)=f(x_1)(t), where t corresponds to the time point.

We denote the Boolean function f with
``
tree<-list(1,c(2,8),c(4,5))
``
The integer 1/2/4 smaller than 5 (the number of nodes) represent the 1st/2nd/5th node respectively; and the integer 8 greater than 5 represents the (8-5)rd node. 

We present the prime implicant using the same way
 ```
 # A prime implicant denotes x1&!x2
 pi<-list(c(1,7))
 ```

<b>Network inference</b>:
- Step. 1 Generate a random Boolean network.
  We use the R package <a href="https://cran.r-project.org/web/packages/BoolNet/index.html">BoolNet</a> to generate a random Boolean network.
  ```
  # load BoolNet and ATEN to library
  library("BoolNet")
  library("ATEN")
  
  # ngenes is the number of nodes in a Boolean network
  ngenes<-10
  # k is the maximum number of genes 
  k<-5
  #call generateRandomNKNetwork to generate a Boolean network, for more deatils about arguments, see BoolNet
  set.seed(0)
  net1<-generateRandomNKNetwork(ngenes, k, topology="scale_free",simplify=TRUE,readableFunctions=TRUE)
  
  # see net1
  > net1
  Boolean network with 10 genes
  
  Involved genes:
  Gene1 Gene2 Gene3 Gene4 Gene5 Gene6 Gene7 Gene8 Gene9 Gene10

  Transition functions:
  Gene1 = (!Gene1)
  Gene2 = (Gene8)
  Gene3 = (Gene10)
  Gene4 = (Gene10)
  Gene5 = (Gene2 & Gene3)
  Gene6 = (Gene4)
  Gene7 = (!Gene2) | (!Gene5)
  Gene8 = (!Gene8) | (Gene4)
  Gene9 = (Gene5)
  Gene10 = (!Gene3)
  ```
- Step. 2 Build the time-series data; the datalist is saved in a list as well.
  ```
  # The time-series data will be generated based on the Boolean functions in the network 'net1'
  # numSeries and numPoints refer to the the number of time series and the number of time points, respectively.
  # noiseLevel represents the that a gene state can randomly flipping with the probability (noiseLevel*100)%
  datalist<-buildTimeSeries(network=net1,numSeries=10,numPoints=10,noiseLevel=0)
  
  ```
- Step. 3 Select a target node, generate the bootstrap samples and out-of-bag (oob) samples for inferring and selecting prime implicants (PIs)
  ```
  # For instance, we select the 3rd node as the target node
  target<-6
  
  # Generate the bootstrap samples and oob samples according to the time-series data
  datasamples<-bootstrap(datalist)
  # respinbag and respoutbag save the expression values of the target node
  datasamples$respinbag<-matrix(datasamples$respinbag[,target])
  datasamples$respoutbag<-matrix(datasamples$respoutbag[,target])
  ```
- Step.4 Initilize the parameters used in ATEN
  ```
  # startT, endT and maxIter refer to the upper, lower temperature (on a log10 scale) and the maximum number of iterations used in simulated annealing algorithm
  # maxK represents the maximum number of input nodes of the target node (the maximum number of leaves in an And/Or tree), it is required when prior knowledge, e.g. the in-degree is 8, is available. If such information is not known, then it can be set as a very large value, e.g. '.Machine$integer.max'
  # rate represents how many non-important PIs are removed in each recursion.
  # nodes represents the number of node in the Boolean network
  parameters<-c(startT=2,endT=-2,maxIter=5000,maxK=5,rate=0.2,nodes=ngenes)
  
  # We shall discuss how to tune those arguments later.
  ```
- Step. 5 Find the important PIs 
  ```
  # load parallel computation to library
  library("parallel")
  
  # B represents the how many trees would be generated
  # pars is the argument for parallel computation
  # the relevant datalist and datasamples are also required for network inference
  # the last parameter 'seed' of PIs() is for helping reproduce the results, we set it as 0 here.
  PIs<-findPIs(B=10,datalist,datasamples,parameters,0)
  
  
  # In our case, we obtained 5 prime implicants after removing non-important ones
  > PIs
  [[1]]
  [1] 4

  [[2]]
  [1] 3

  [[3]]
  [1]  4 19

  [[4]]
  [1]  4  8 19

  [[5]]
  [1]  4 13 18

  ```
  -Step. 6 Find the Boolean function according to those PIs and RFRE framework
  ```
  # Update the datasets so that the new datalist corresponds to the PIs
  datalist[[2]]<-generateData(PIs,datalist)
  datalist[[3]]<-matrix(datalist[[3]][,target])
  datasamples<-bootstrap(datalist)
  datasamples$respinbag<-matrix(datasamples$respinbag)
  datasamples$respoutbag<-matrix(datasamples$respoutbag)
  
  # Identify the final Boolean function
  # the last parameter of findBF() is for helping reproduce the results, we set it as 0 here.
  BF<-findBF(10,PIs,target,parameters,datalist,datasamples,0)
  
  # Check the final solution we obtained
  >BF
  "Gene4"
  ```
  <b>Something new in ATEN</b>. In some cases, Step. 6 is not required, for instance, using the same Boolean network but with noisy data this time
  ```
  # Generate the time-series data with 5% noise
    datalist<-buildTimeSeries(network=net1,numSeries=10,numPoints=10,noiseLevel=0.05)
  
  # Now selected the fifth node as the target node
  target<-1
  
  # Generate the bootstrap samples and oob samples according to the time-series data
  datasamples<-bootstrap(datalist)
  # respinbag and respoutbag save the expression values of the target node
  datasamples$respinbag<-matrix(datasamples$respinbag[,target])
  datasamples$respoutbag<-matrix(datasamples$respoutbag[,target])
  
  # Find the important PIs
  PIs<-findPIs(B=10,datalist,datasamples,parameters,0)
  
  # See PIs
  > PIs
  [1] "!Gene1"
  
  # We can find the result is not a list of PIs but the final Boolean function. 
  # And the final Boolean function  "Gene10&!Gene4" contains the true input node (i.e. Gene10) and a false input node (i.e. Gene4). 
  # In this way, Step.6 is not required.
  # There are two reasons:
  # 1) only 1 PI left after eliminating the non-important ones;
  # 2) after elimination, only a few PIs (=<4) left, we then directly build a And/Or tree.
  # Think about the 2nd issue, actually a better way is to invoke Best-Fit method there to find the optimal solution. 
  # And also, we can enlarge the threshold (e.g. 8) as Best-Fit can always handles such small networks with good performance.
  ```

<b>Parameters setting</b><br/>
  In our opinion, the default parameters/arguments values are OK for small networks (<=10 nodes), but not the best; and you can make them better if you're willing to invest time in learning how to set the parameters/arguments. 
  - B, the number of And/Or trees in ATEN. This parameter always depends on the size/noise of the data. We suggest you set a higher B for large/noisy networks. If you set a very large B, although it won't affect the results much, but it takes much time... 
  - startT, the uppper temperature is not suggested to be very high, as more false pis would be introduced; and if it is too low, then only a 'local' solution space will be searched. But a higher termperature is still better than a lower one. Normally, we set it as 2.
  - endT, as the probability of acceptation depends on the temperature, we suggest tune it according to the number of acception iterations. You can print the acceptances/the quality of solutions if you want. In general, if there are many acceptances, lower end is better; if no acceptance in a row, raise it a bit.  Normally, we set it as a negative value less than -2. 
   Note that we run tree growth at fixed temperature until the iteration exceed a threshold. It also can be changed. Users can try to use another cooling scheme. One idea is to adjust the cooling scheme to the number of acceptances. And also any other good SA algorithm is also welcome to be introduced (which we will make it soon).
  - maxIter, the number of iterations. This iterations actually is related to the temperature and the scale of the dataset. If only 5 input nodes are included, then 10000 iterations (>>2^2^5) are not reasonable and therefore we are thinking about introducing Best-fit method. Normally we suggest 10000 iterations if the network size is greater than 10 (2^2^10 possible Boolean functions). And if you find there are still many good solutions in the end, then you may need to increase the number of iterations and tune the endT as well. By default we set it as a value between 5000 and 10000.
  By the way, it is very easy to make ATEN as a feature selection tool before applying Best-Fit (i.e. finding all putative Boolean functions). We shall update it later.
  - However, we have to admit that this investment takes much time... Different sizes of networks have differnt in-degree. And also the number of PIs always changes in the RFRE framework, and it is not possible to tune parameters for each individual node. Hence in a word, we suggest using a higher startT in the beginning, and probably reduce it and the maxIter according to the result returned of each recursion. We did not try to identify what a thresold is the best, but we think 10000 iterations for >10 nodes performs good.  
  
 Someone would be interested in how to set the tree size. We have discuss more details in our Supplementary Data. In a word, the true positive rate (TPR or Recall) would increase gradually and become stable when the tree size is greater than a specific threshold. And the false postive rate would increase gradually as more and more false positive interactions would be included.
  
  
 <b>Future work</b><br/> 
  - Besides what we discuss above, another direction is to make it adaptive for different sizes of networks. 
  - We also expected our idea can be used for inferrng probabilistic Boolean networks and asynchronous networks.
  
<br/>


### References: <br/>
---
*Shi, N., Zhu. Z.X., Tang, K., Parker, D. and He, S., 2019. And/Or Tree Ensemble for inferring accurate Boolean network topology and dynamics*

*Müssel, C., Hopfensitz, M. and Kestler, H.A., 2010. BoolNet—an R package for generation, reconstruction and analysis of Boolean networks. Bioinformatics, 26(10), pp.1378-1380.*

*Lähdesmäki, H., Shmulevich, I. and Yli-Harja, O., 2003. On learning gene regulatory networks under the Boolean network model. Machine learning, 52(1-2), pp.147-167.*
