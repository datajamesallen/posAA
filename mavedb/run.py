#!/usr/bin/python3

"""
main execution entry, used to manage the project
"""

def default_mode():
    print('# Updating analysis')
    import analysis
    print('# Adding extra information')
    import add_info
    print('# Performing meta analysis')
    import meta
    print('# Done')
    return None

default_mode()
