#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
USPEX 9.4.4 release
2015 Oganov's Lab. All rights reserved.
"""

import os


# -------------------------------------------------------------------------------
def uspex_vars(envvar):
    if envvar in os.environ:   # modified on 210421
#    if os.environ.has_key(envvar):
        return_value = None
    else:
        return_value = envvar + ' environment variable is not set. Exit.'

    return return_value
