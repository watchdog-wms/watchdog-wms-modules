'''
Created on 27 Apr 2017

@author: friedl
'''


class Table(object):
    '''
    representation of a table
    '''
    
    def __init__(self, reftable=None):
        '''
        Constructor
        'reftable' is another instance of Table, if it is not None, the constructor copies the column specifications of reftable
        '''
        
        # table size
        self._rows = 0
        self._columns = 0
        
        # column specifications: column names and stored data types
        self._column_to_Pos = {}
        self._column_names = []
        self._types = []
        # copy column specifications from reftable
        if(reftable is not None):
            self._columns=reftable._columns
            self._column_to_Pos = reftable._column_to_Pos.copy()
            self._column_names = reftable._column_names.copy()
            self._types = reftable._types.copy()
        
        # list of row lists
        self._content = []
        
    def copy(self):
        '''
        creates a copy of the current table
        '''
        
        copied = Table(self)
        for row in range(0, self.rowNum()):
            row_copy = list(self.getRow(row))
            copied.addRow(row_copy)
        return copied
        
    def size(self):
        '''
        returns the size of the table
        '''
        
        return(self._rows, self._columns)
    
    def rowNum(self):
        '''
        returns the number of rows in the table
        '''
        
        return self._rows
    
    def colNum(self):
        '''
        returns the number of columns in the table
        '''
        
        return self._columns
    
    
    def addColumn(self, columnType, columnName=None, defaultValue=None):
        '''
        adds a column at the end of the table, extends all rows with defaultValue
        '''
        
        # adapt column data structures
        self._types.append(columnType)
        if(columnName is not None):
            if(columnName not in self._column_names):
                self._column_to_Pos[columnName]=self._columns
            # prevent duplicate column names
            else:
                raise ValueError('Cannot add column '+columnName+': duplicate column name!')
        self._column_names.append(columnName)
        self._columns+=1
        
        # extend rows with default value
        for row in self._content:
            row.append(defaultValue)
            
        
    def addRow(self, rowcontent):
        '''
        adds a row at the end of the table, content is a list of values of any type
        '''
        
        if(len(rowcontent)==self._columns):
            rowToAdd = []
            for cell, typeFunc in zip(rowcontent, self._types):
                rowToAdd.append(typeFunc(cell))
            self._content.append(rowToAdd)
            self._rows+=1
        else:
            raise ValueError('Row length mismatch')
        
    
    def changeColumnType(self, column, coltype):
        '''
        changes the type of a column
        '''
        
        colPos = self._check_table_access_column(column, 'changeColumnType')
        self._types[colPos] = coltype
        for row in self._content:
            row[colPos] = coltype(row[colPos])
    
            
    def changeColumnName(self, oldName, newName):
        '''
        changes the name of a column
        'oldName' gives the column to rename as position or name
        '''
        
        # check new name
        if(newName in self._column_to_Pos):
            raise ValueError('Cannot rename column '+oldName+' to '+newName+': Duplicate column name!')
        
        # check column to rename
        colPos = self._check_table_access_column(oldName, 'ChangeColumnName')
        
        # perform rename: update column -> pos dictionary and column name list
        self._column_to_Pos[newName]=colPos
        nameToRemove = self._column_names[colPos]
        self._column_names[colPos]=newName
        if(nameToRemove is not None):
            del self._column_to_Pos[nameToRemove]


    def getColumnNames(self, noneHandling=True):
        '''
        returns a list of the column names, the position of the name within the list correspond to the position of the column
        noneHandling specifies handling of nameless column:
        true <-> columns without name get the NLC[1-9]+ (nameless column no. [1-9]+)
        false <-> return None as name
        '''
        
        colNames=[]
        none_counter=1
        for col_name in self._column_names:
            if(col_name is None and noneHandling is True):
                colNames.append('NLC'+str(none_counter))
                none_counter+=1
            else:
                colNames.append(col_name)
        
        return colNames
    
    def getColumnName(self, column):
        '''
        returns the name of the column 'column'
        ATTENTION: in contrast to getColumnNames this function returns always None if the column does not have name
        '''
        
        colPos = self._check_table_access_column(column, 'GetColumnType')
        return self._column_names[colPos]
    
    
    def getColumnType(self, column):
        '''
        returns the type of the column 'column'
        '''            
        
        colPos = self._check_table_access_column(column, 'GetColumnType')
        return self._types[colPos]
    
    
    def get(self, row, column):
        '''
        get a value of the table
        '''
        
        rowPos = self._check_table_access_row(row, 'Get')
        colPos = self._check_table_access_column(column, 'Get')
        
        return self._content[rowPos][colPos]


    def getColumn(self, column):
        '''
        returns a list with the values in the column, the position within the returned list corresponds to the row number of the value
        '''
        
        colPos = self._check_table_access_column(column, 'GetColumn')
        
        res = []
        for row in self._content:
            res.append(row[colPos])
        return res
    
    
    def getRow(self, row, select_cols=None):
        '''
        returns a list (copy) with the values in the row, the position within the returned list corresponds to the column number of the value
        '''
        
        # create a copy of the row as return
        res = []
        rowPos = self._check_table_access_row(row, 'GetRow')
        
        # return complete row
        if(select_cols is None):
            for col in self._content[rowPos]:
                res.append(col)
        
        # return specific columns of the row
        else:
            rowContent = self._content[rowPos]
            for resCol in select_cols:
                colPos = self._check_table_access_column(resCol, 'GetRow')
                res.append(rowContent[colPos])
        
        return res
    
    
    def set(self, row, column, value):
        '''
        set a value of the table
        '''
        
        rowPos = self._check_table_access_row(row, 'Set')
        colPos = self._check_table_access_column(column, 'Set')
        
        # check data type of column
        if(not isinstance(value, self._types[colPos])):
            value = self._types[colPos](value)
        
        self._content[rowPos][colPos] =value
        
    
    def modifyColumn(self, column, modifying_function, wholeRow=False):
        '''
        applies a function to the value of 'column' in all rows
        modifying_function: old value -> new value
        wholeRow = True -> modifying_function takes as input the table and the current row index
        '''
        
        colPos = self._check_table_access_column(column, 'ModifyColumn')
        
        # modify value of the column based on the old value of the column
        if wholeRow is False:
            for row in self._content:            
                row[colPos] = modifying_function(row[colPos])
        
        # modify value of the column based on the whole current row
        else:
            for rowPos, row in enumerate(self._content):
                row[colPos] = modifying_function(self, rowPos)
            
    #TODO: method to modify a column based on the values of several other columns -> extend method above with list of columns to use
    
    #TODO: WARNING not tested yet completely!
    def sortRows(self, columns, columnKeys, ascending):
        '''
        method to sort the table IN PLACE!
        'columns': list of columns used for sorting, order of column names or positions specifies sorting priority
        'columnKeys': list of functions colum value -> key used for sorting
        'ascending': list of booleans, True -> ascending order, False -> descending order
        '''
        
        for column, keyFunc, asc in reversed(list(zip(columns,columnKeys, ascending))):
            colPos = self._check_table_access_column(column, 'SortRows')
            sortFunc = lambda row: keyFunc(row[colPos])
            if asc is True:
                self._content.sort(key=sortFunc, reverse=False)
            else:
                self._content.sort(key=sortFunc, reverse=True)
            

    
    def _check_table_access_row(self, row, operationString):
        '''
        helper method for getter and setter methods using a row coordinate
        '''
        
        # check row coordinate
        if((not isinstance(row, int)) or row>=self._rows):
            raise ValueError(operationString+'outside table dimensions: row '+repr(row))
        
        return row
    
    
    def _check_table_access_column(self, column, operationString):
        '''
        helper method for getter and setter methods using a column coordinate
        '''
        
        # convert column name to column pos
        colPos=column
        if(isinstance(column, str)):
            if(column in self._column_to_Pos):
                colPos=self._column_to_Pos[column]
            else:
                raise ValueError(operationString+'unknown column name: '+repr(column))
        
        # check column coordinate
        if((not isinstance(colPos, int)) or colPos>=self._columns):
            raise ValueError(operationString+'outside table dimensions: column '+repr(colPos))
        
        return colPos
    