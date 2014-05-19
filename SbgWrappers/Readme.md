This directory contains wrapper code that allows Mitty components to run on the Seven Bridges Platform. Instructions
for wrapping are found here and might be helpful to devs.

Setting up, directory structure etc. etc.
----------
_I don't understand docker/VM/vagrant internals. My user's view of what we are doing as we write wrappers for
our tools and then upload these wrappers to the Seven Bridges computing platform is set out as annotations in
this recipe book. They will probably be wrong, please tell me so I can correct the statements._

Please read the initial part of the [SBG pipeline documentation][sbgreadme] to see how to install Vagrant and
VirtualBox.

[sbgreadme]: http://pythonhosted.org/sbgsdk/tutorial.html


### Setting up the SBG virtual machine.
Seven Bridges supplies us with a virtual machine image (an machine running Ubuntu linux). The supplied archive looks
something like this

    .
    |____.DS_Store
    |____sbg_wrappers
    | |____.DS_Store
    | |____.sbdk
    | | |____logging.json
    | | |____state.json
    | |____head.json
    | |____sbg_wrappers
    | | |______init__.py
    | | |____head.py
    | | |____samtools_view.py
    | |____setup.py
    | |____test-data
    | | |____file.txt
    | | |____mock.sam
    |____Vagrantfile

The tutorial tells us to


Vagrant is a useful
tool for configuring virtual machines such that we can go into the command line and then type `vagrant ssh` to
pretend we are `ssh`-ing into a virtual machine, when in actuality we are starting up an instance of the VM
and



In my case, where both the tool code and the wrappers are under development I like to

  a. Be able to do a `git pull` to grab the latest tool code
  b. Be able to edit the wrapper code and have it reflected in the VM, and
  c. Have the wrappers in the same repository as the tool, under a separate directory

-->>> note on where to place .sbdk and Vagrant for this

Place the downloaded Vagrant image file where you wish (In my case I have it at `~/Platform` which I changed from the
default `sbg-projects`). There should be a Vagrant image vile (called Vagrant).

Now start the virtual machine with the VirtualBox image with the Seven Bridges SDK from within this directory.

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

    sbg init mitty  # There might be a noticable pause here. You will be asked for your SBG Platform credentials
                    # Nebojsa says sbg is going to the platform to make sure there are no name clashes

Choose this name well because this is what you will see on the Apps panel.

We need to be in this directory to run any `sbg` commands.

    cd mitty        #
    sbg sh          # enter a docker container where we can install Mitty

Now, install the dependencies that Mitty needs that is not included with the sbg image. This will take a while.

    pip install docopts pysam pyvcf numpy

Now get out of the container to commit

    exit

**Interestingly, any `sbg` command other than `sbg sh` seems to require an internet connection. I'm not sure why. This
puts a mild damper on things, since I often write code on the train etc. where I don't have an internet connection.**

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

    sbg test mitty.plugins.mutation.snp_wrapper
    sbg test mitty.plugins.mutation.insert_wrapper

An interesting command to run from within a container is

    sbg schema

This reads the wrappers listed in `__init__.py` and drops json output that looks like:

      {
        "schema": {
          "inputs": [],
          "outputs": [
            {
              "description": "Deletion definitions for mutate",
              "id": "json_fragment",
              "list": false,
              "name": "Deletions",
              "required": false,
              "types": []
            }
          ],
          "params": [
            {
              "category": "General",
              "condition": null,
              "default": null,
              "description": "A unique name for this instance of the insert generator",
              "id": "model_id",
              "list": false,
              "name": "Model id",
              "pattern": null,

                  ...

              "required": false,
              "step": null,
              "type": "integer"
            }
          ]
        },
        "wrapper_id": "mitty.plugins.mutation.delete_wrapper.Deletion"
      }


Push the wrappers

    sbg push "My commit message which will become the name for this release of the Mitty toolkit"

I usually end up doing something like

    sbg push "v0.1.2"

because that's what appears in the Apps panel on the left and I want that to be informative.

If you get `No wrappers registered (empty __init__.py?). Exiting.` then you need to go fill out `__init__.py` at the
root level. For example

    from .mutate_wrapper import Mutate
    from .plugins.mutation.snp_wrapper import SNP
    from .plugins.mutation.insert_wrapper import Insert
    from .plugins.mutation.delete_wrapper import Deletion

This is cool because it allows you to select which wrappers get pushed.

Check at this url to see if the image took

    https://images.sbgenomics.com/v1/repositories/kghosesbg/mitty/images



