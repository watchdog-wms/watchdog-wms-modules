### Contributing

Contributions to the [watchdog-wms-modules](watchdog-wms/watchdog-wms-modules) repository are welcome!
If you want to share a module, please follow these steps:

1) fork the repository and create a new branch for your new module
2) commit your module into your fork
    - either use the webinterface of Github (buttons are located left of green "Clone or Download" button)
        - Create new file: create and commit a single file in a web editor
        - Upload Files Button: use drag & drop gestures to upload multiple files or complete folders
    - or clone the repository and use the git command line tools to create and push a new commit
3) create a pull request
4) automatic tests will be performed on your pull request (by Travis CI)
5) members of the contributor team can merge your pull request if all checks succeed

More detailed instructions how to contribute to a github repository can be found [here](https://github.com/firstcontributions/first-contributions).


### Join the contributor team
Members of the contributor team have write access to the repository and can merge pull requests if all Travis CI checks succeeded.

To become a member, please follow these steps:
1) share a few modules
2) request to become a member by mentioning watchdog-wms-bot in one of your pull requests (just include '@watchdog-wms-bot' in your comment)

### Module structure

Modules must contain at least an XSD module definition file (*moduleName.xsd*) and an XML documentation file (*moduleName.docu.xml*). 
If test data is included, it should be stored in a folder named *test_data*.
Files shared between modules must be located in the [sharedUtils/](https://github.com/watchdog-wms/watchdog-wms-modules/tree/master/sharedUtils) folder.

Example:

    moduleName
    ├── moduleName.xsd
    ├── moduleName.docu.xml
    ├── moduleName.sh
    └── test_data
        ├── testFile1.txt
        └── testFile2.txt
    
More information on how to create a module can be found in the [documentation](https://rawgit.com/klugem/watchdog/master/documentation/Watchdog-manual.html#custom_modules).
   
### Automatic tests on pull requests

Before a pull request can be merged with the master branch some Travis CI tests must succeed.
Currently the following tests are implemented:

- signed commit test: 
  - only [signed commits](https://help.github.com/en/articles/signing-commits) can be merged 
- write permission test:
  - the github username of the creator of the pull request must be included in the XML documentation file (*\<maintainer\>* tag)
- separate module test: 
  - all files affected by the pull request must be part of one module (exception: files located in [sharedUtils/](https://github.com/watchdog-wms/watchdog-wms-modules/tree/master/sharedUtils) are allowed)
- XSD validation test: 
  - XSD module definition file must be parseable by Watchdog
- XML documentation test:   
  - XML documentation file must follow its [XSD schema](https://github.com/watchdog-wms/watchdog-wms/blob/master/xsd/documentation.xsd)
  - parameter and return values must fit 
- virus scanner test: 
  - all files part of the pull request must pass the virus scan
