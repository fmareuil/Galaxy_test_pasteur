from module_env import *


modulesliste = moduleavail()
for module in modulesliste:
    if "(default)" in module:
        indexe = modulesliste.index(module)
        modulesliste.pop(indexe)
        modulesliste.insert(indexe,module.replace("(default)",""))


dictmodule = {}
for module in modulesliste:
    commandsmodule = modulehelp(module).split()
    try:
        commandsmodule = commandsmodule[commandsmodule.index("commands:")+1:]
    except:
        print module
    dictmodule[module] = commandsmodule
    
command = "fastx_collapser"
for module in dictmodule:
    if command in dictmodule[module]:
        print module
 
