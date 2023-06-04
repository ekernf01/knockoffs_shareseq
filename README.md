## SHARE-seq skin experiments

Code to repeat our experiments on the SHARE-seq mouse skin data for our manuscript "Model-X knockoffs reveal data-dependent limits on regulatory network identification". For more info, see the [central project repo](https://github.com/ekernf01/knockoffs_paper). Anyone should be able to repeat our experiments by using the same environment we used on AWS (Ubuntu 18.04); see `run_on_aws.sh` for details. Plots should be saved to a sub-directory called `v12`.

If you have trouble setting up the environment, you can also use our docker image:

```sh
docker pull    ekernf01/knockoffs_shareseq
docker run -it ekernf01/knockoffs_shareseq
```

You probably will need to run those commands using sudo. You may need to make a couple more changes too: in our limited tests, most r packages were installed but ours, rlookc, was merely downloaded and not installed. That's fixable with `install.packages("rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))`, and we just don't want to make this change in the install script because it might break backwards compatibility with our current repro instructions. Docker is not how we ran our experiments initially, but with AWS already phasing out Ubuntu 18.04, we hope this Docker image will allow exact repro further into the future. 
