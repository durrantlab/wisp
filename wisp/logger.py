def log(astring, fileobjects):  # prints to screen and to log file
    """Outputs WISP messages

    Arguments:
    astring -- a string containing the message
    fileobjects -- a list of python file objects specifying where the messages should be saved
    """

    if not isinstance(fileobjects, list):
        # it's not a list, so make it one
        fileobjects = [fileobjects]

    print(astring)

    for fileobject in fileobjects:
        fileobject.write(astring + "\n")
