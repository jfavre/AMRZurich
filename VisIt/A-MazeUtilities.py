# version 1.00
# written by Jean M. Favre, CSCS
# tested with VisIt 2.4.2
# last change: Sun Mar 18 12:33:40 CET 2012
########################################################################
# The AMAZE reader for VisIt gives names for levels and patches
# Levels in the AMR hierarchy start at 0
# Patches start at 1
#
# When a Plot is selected, it is "Active". The subset menu in the GUI lets
# you turn ON/OFF different levels and patches
#
# The following Python functions do the same for Levels or patches
# Example:
#
# 3 Levels names "level0", "level1", "level2"
#
# Level0 has one patch called "Grid 1"
# Level1 has two patches called "Grid 2" and "Grid 3"
# Level2 has two patches called "Grid 4" and "Grid 5"
#
# type "TurnOffLevel(0)" or "TurnOffPatch(2)"
#
# Additionally, the function ListStars() will give names and indices of stars,
# such that, giving an index, a star can be turned OFF/ON with codes similar:
# silr = SILRestriction()
# silr.TurnOnSet(120)
# SetPlotSILRestriction(silr ,0)
# The same is done for star names, e.g.
#
# ListStars()
# TurnOnStar('IsotrInfWind')
# TurnOnStar('SimpleAccretor')
#
# TurnOffStar('IsotrInfWind')
# TurnOffStar('SimpleAccretor')
#
########################################################################

def TurnOffLevel(level):
   levname = "level%d" % (level)
   # The levels category has the substring "Level" in it
   s = SILRestriction()
   cats = s.Categories()
   level_category = ""
   for x in cats:
      if "Level" in x:
          level_category = x
   if level_category == "":
       print "*** Unable to locate level category ***\n"
       return
   levels = s.SetsInCategory(level_category)
# Visit numbers "level sets" starting at 1 (which is our "level0") in the list
   foundOne = False
   for l in levels:
       levnum = s.SetName(l)
       if levnum == levname:
           s.TurnOffSet(l)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate level %d ***\n" %(level)
       return
   SetPlotSILRestriction(s)

def TurnOnLevel(level):
   levname = "level%d" % (level)
   # The levels category has the substring "Level" in it
   s = SILRestriction()
   cats = s.Categories()
   level_category = ""
   for x in cats:
      if "Level" in x:
          level_category = x
   if level_category == "":
       print "*** Unable to locate level category ***\n"
       return
   levels = s.SetsInCategory(level_category)
# Visit numbers "level sets" starting at 1 (which is our "level0") in the list
   foundOne = False
   for l in levels:
       levnum = s.SetName(l)
       if levnum == levname:
           s.TurnOnSet(l)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate level %d ***\n" %(level)
       return
   SetPlotSILRestriction(s)

def TurnOffPatch(patch):
   gridname = "Grid %d" % (patch)
   # The patch category has the substring "Patches" in it
   s = SILRestriction()
   cats = s.Categories()
   patch_category = ""
   for x in cats:
      if "Patches" in x:
          patch_category = x
   if patch_category == "":
       print "*** Unable to locate patch category ***\n"
       return
   patches = s.SetsInCategory(patch_category)
# Visit numbers "patch sets" starting at 1 (which is our "Grid 1") in the list
   foundOne = False
   for p in patches:
       patnum = s.SetName(p)
       if patnum == gridname:
           s.TurnOffSet(p)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate patch %d ***\n" %(patch)
       return
   SetPlotSILRestriction(s)

def TurnOnPatch(patch):
   gridname = "Grid %d" % (patch)
   # The patch category has the substring "Patches" in it
   s = SILRestriction()
   cats = s.Categories()
   patch_category = ""
   for x in cats:
      if "Patches" in x:
          patch_category = x
   if patch_category == "":
       print "*** Unable to locate patch category ***\n"
       return
   patches = s.SetsInCategory(patch_category)
# Visit numbers "patch sets" starting at 1 (which is our "Grid 1") in the list
   foundOne = False
   for p in patches:
       patnum = s.SetName(p)
       if patnum == gridname:
           s.TurnOnSet(p)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate patch %d ***\n" %(patch)
       return
   SetPlotSILRestriction(s)

def ListStars():
  for s in SILRestriction().SetsInCategory('domains'):
    print s, ": ", SILRestriction().SetName(s)

def TurnOffStar(starname):
   # The category should be "domains"
   s = SILRestriction()
   cats = s.Categories()
   dom_category = ""
   for x in cats:
      if "domains" in x:
          dom_category = x
   if dom_category == "":
       print "*** Unable to locate 'domains' category ***\n"
       return
   stars = s.SetsInCategory(dom_category)
# Visit gives names to the stars as found in the HDF5 file
   foundOne = False
   for star in stars:
       mystar = s.SetName(star)
       if mystar == starname:
           s.TurnOffSet(star)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate star ", starname
       return
   SetPlotSILRestriction(s)

def TurnOnStar(starname):
   # The category should be "domains"
   s = SILRestriction()
   cats = s.Categories()
   dom_category = ""
   for x in cats:
      if "domains" in x:
          dom_category = x
   if dom_category == "":
       print "*** Unable to locate 'domains' category ***\n"
       return
   stars = s.SetsInCategory(dom_category)
# Visit gives names to the stars as found in the HDF5 file
   foundOne = False
   for star in stars:
       mystar = s.SetName(star)
       if mystar == starname:
           s.TurnOnSet(star)
           foundOne = True
   if not foundOne:
       print "*** Unable to locate star ", starname
       return
   SetPlotSILRestriction(s)

 
