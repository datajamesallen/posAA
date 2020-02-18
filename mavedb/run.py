#!/usr/bin/python3

"""
main execution entry, used to manage the project
"""

def default_mode():
    print('# Updating analysis')
    import analysis
    print('# Performing meta analysis')
    import meta
    print('# Done')
    return None

default_mode()
