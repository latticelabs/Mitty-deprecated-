This directory contains wrapper code that allows Mitty components to run on the Seven Bridges Platform. Instructions
for wrapping are found here and might be helpful to devs.

Setting up
----------
Please read the initial part of the [SBG pipeline documentation][sbgreadme] to see how to install Vagrant and
VirtualBox.

[sbgreadme]: http://pythonhosted.org/sbgsdk/tutorial.html

Place the downloaded Vagrant image file where you wish (In my case I have it at `~/Platform` which I changed from the
default `sbg-projects`). There should be a Vagrant image vile (called Vagrant).


Download the VirtualBox image with the Seven Bridges SDK from within this directory.

    vagrant up

This creates a `.vagrant` directory in which the virtual machine state is housed.

Also note that the log information that Vagrant puts out indicates that (for example in my case) the virtual machine
folders are mapped to the outside.

    /vagrant => /Users/kghose/Platform
    /home/vagrant/projects => /Users/kghose/Platform

Log in to this virtual machine and pull the docker image.

    vagrant ssh                         # Start an ssh session with the guest machine.
    docker pull sevenbridges/sdkbase:beta  # I don't know why I need to do this, but this is the magic sauce to get us started.

Now, update the SBG SDK

    sudo pip install sbgsdk --upgrade   # Make sure we have the latest SDK.


Now, create and enter a docker container where we can install Mitty and its dependencies. If you place this under
`/vagrant` or `/home/vagrant/projects` you will be able to see the directory from your (regular) host machine. This
has great advantages when developing the code for the wrappers.

    cd /home/vagrant/projects

    sbg init Mitty  # There might be a noticable pause here. You will be asked for your SBG Platform credentials
                    # Nebojsa says sbg is going to the platform to make sure there are no name clashes

We need to be in this directory to run any `sbg` commands.

    cd Mitty        #
    sbg sh          # enter a docker container where we can install Mitty

Now, install the dependencies that Mitty needs that is not included with the sbg image. This will take a while.

    pip install docopts pysam pyvcf numpy

Now get out of the container to commit

    exit

If you execute the `docker images` command you will get a list of docker images. The one you just committed will show up
tagless and repositoryless, something like

    REPOSITORY                              TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
    <none>                                  <none>              8ab3de688526        32 seconds ago      614.7 MB
    images.sbgenomics.com/kghosesbg/mitty   a7ee5ec45d33a745    a7ee5ec45d33        59 minutes ago      560.2 MB
    sevenbridges/sdkbase                    beta                60c297c6f611        11 weeks ago        560.2 MB


Not needed but we can ...
... I like to tag the version with all the tools installed in case I need to revert to it for some reason.

    docker tag 8ab3de688526 images.sbgenomics.com/kghosesbg/mitty:base

Now `docker images` will look something like

    REPOSITORY                              TAG                 IMAGE ID            CREATED              VIRTUAL SIZE
    images.sbgenomics.com/kghosesbg/mitty   base                8ab3de688526        About a minute ago   614.7 MB
    images.sbgenomics.com/kghosesbg/mitty   a7ee5ec45d33a745    a7ee5ec45d33        About an hour ago    560.2 MB
    sevenbridges/sdkbase                    beta                60c297c6f611        11 weeks ago         560.2 MB

Now check out the latest version of Mitty. We should be sitting under

    rm -rf mitty/   # I like to use git to download Mitty and sbg init pre creates a folder
    git clone https://github.com/latticelabs/Mitty.git  # You will be asked for login and password since this is a closed repo.


    sbg test snp_wrapper


    sbg push "My commit message which will become the name for this release of the Mitty toolkit"


Check at this url to see if the image took

    https://images.sbgenomics.com/v1/repositories/kghosesbg/mitty/images

