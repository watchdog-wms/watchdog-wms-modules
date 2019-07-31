### Folder structure

Modules must contain at least an XSD module definition file (*moduleName.xsd*) and an XML documentation file (*moduleName.docu.xml*).
If test data is included, it should be stored in a folder named *test_data*.

Example:

    moduleName
    ├── moduleName.xsd
    ├── moduleName.docu.xml
    ├── moduleName.sh
    └── test_data
        ├── testFile1.txt
        └── testFile2.txt
   
### Automatic Travis CI tests triggered by pull requests

Before a pull request can be merged with the master branch some Travis CI tests must succeed.
Currently the following tests are implemented:

- signed commit test: 
  - only [signed commits](https://help.github.com/en/articles/signing-commits) can be merged 
- write permission test:
  - the github username of the creator of the pull request must be included in the XML documentation file (*\<maintainer\>* tag)
- separate module test: 
  - all files affected by the pull request must be part of one module
- XSD validation test: 
  - XSD module definition file must be parseable by Watchdog
- XML documentation test:   
  - XML documentation file must follow its [XSD schema](https://github.com/watchdog-wms/watchdog-wms/blob/master/xsd/documentation.xsd)
  - parameter and return values must fit 
- virus scanner test: 
  - all files part of the pull request must pass the virus scan
