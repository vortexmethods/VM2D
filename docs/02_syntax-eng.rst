Source files syntax 
===================

All files with source data and settings for **VM2D** are text files and have the uniform *C++*-like syntax.

*      Two types of comments are allowed, following the *C++* pattern:

       *     ``//`` - comment from current position to end of line;

       *     ``/* ... */`` - comment of a fragment of arbitrary length, including multi-line (such comments cannot be nested);

*      line break is equivalent to a space;
*      spaces in source files and tabs are ignored. 


The source data files are organized as *dictionaries*, i.e. according to principle

.. code-block:: c
   :name: key1
   
   key1 = value1;
   key2 = value2;
   ...

*Keys* are case insensitive, i.e. ``key1``, ``Key1``, ``KEY1``, ``kEy1`` are equal.

*Values* can be either *simple* values or *lists*, if the latter, list is enclosed in curly braces, and values in the list are separated by commas

.. code-block:: c
   :name: key2
   
   key1 = value1;
   key2 = {value21, value22};
   key3 = {value31, value32, value33,};

It is permissible to leave a comma at the end of the list of values before the closing curly brace; it will be ignored in this case.

The type of values is arbitrary --- they can be integers, floating point numbers, logical values (``true``/``false``, or, which is the same, ``yes``/``no`` and ``0``/``1``), symbols, strings. String values without spaces can be specified either in quotes or without them; if the string value contains a space, comma, semicolon or brackets, then it must be enclosed in double quotes:

.. code-block:: c
   :name: key3
   
   integerkey1       = 12;
   floatingpointkey2 = 6.22;
   vectorkey3        = {3.14159265, 2.71828183};
   logicalkey4       = true;
   logicalkey5       = no;
   stringkey6        = "Hello, World!";  // string-type values are
   stringkey7        = Goodbye;          // case sensitive

The type of the values corresponding to the keys is not analyzed or controlled when files are read; the values themselves are initially read simply as sets of characters (strings). When accessing the corresponding parameter from the executable program, the functions for converting the read lines to the required type are called.

Some values can be *complicated*, i.e. they, in turn, depend on some parameters, in this case the parameters of such values are indicated in parentheses

.. code-block:: c
   :name: key4
   
   key1 = value1(parameter1);




