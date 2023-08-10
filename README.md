## SHARE-seq skin experiments

Code to repeat our experiments on the SHARE-seq mouse skin data for our manuscript "General conditional hypothesis tests suggest fundamental limits on regulatory network identification". For more info, see the [central project repo](https://github.com/ekernf01/knockoffs_paper). You can repeat our experiments by using a clean install of Ubuntu 18.04 or Ubuntu 20.04 and running `run_on_aws.sh`. Plots will be saved to a sub-directory called `v12`.

If you have trouble accessing a clean ubuntu install, you can also use our docker image:

```sh
sudo docker pull    ekernf01/knockoffs_shareseq
sudo docker run -it --rm ekernf01/knockoffs_shareseq
```

This will land you in an interactive shell in the Docker container, and you can run our experiments by running `run_in_docker.sh`.
