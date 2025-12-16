vtk-developers [vtk-developers-bounces@vtk.org] on behalf of Shawn Waldon [shawn.waldon@kitware.com]
[Reply] [Reply All] [Forward]
Actions
To:
 Ken Martin ‎[ken.martin@kitware.com]‎ 
Cc:
 vtk-developers ‎[vtk-developers@vtk.org]‎‎; Will Schroeder ‎[will.schroeder@kitware.com]‎ 
 
Friday, January 29, 2016 6:23 PM
Ken,

What you are wanting sounds very similar to how the streaming particles plugin in ParaView works.  In this case it is the representation that is deciding what data to request from the upstream pipeline.  I have an example source in there [1] that randomly generates an octree-like structure that it expects (of vtkMultiBlockDataset).  Either you or Will may want to look at how that is working while you are tackling this.

Shawn

ParaView-src/Plugins/StreamingParticles/vtkPVRandomPointsStreamingSource.{h,cxx}.

