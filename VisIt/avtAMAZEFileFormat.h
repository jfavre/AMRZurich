/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtAMAZEFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_AMAZE_FILE_FORMAT_H
#define AVT_AMAZE_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#define H5_USE_16_API
#include <hdf5.h>

#include <vtkAMAZEReader.h>
#include <void_ref_ptr.h>
#include <vector>

class DBOptionsAttributes;

// ****************************************************************************
//  Class: avtAMAZEFileFormat
//
//  Purpose:
//      Reads in AMAZE files as a plugin to VisIt.
//
//  Programmer: jfavre -- generated by xml2avt
//  Creation:   Fri Nov 6 10:10:36 PDT 2009
//
// ****************************************************************************

class avtAMAZEFileFormat : public avtSTMDFileFormat
{
  public:
                       avtAMAZEFileFormat(const char *, DBOptionsAttributes *);
    virtual           ~avtAMAZEFileFormat() {this->FreeUpResources();};

    virtual const char    *GetType(void)   { return "AMAZE"; };
    virtual void           FreeUpResources(void); 
    virtual void           ActivateTimestep(void);
    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    void *GetAuxiliaryData(const char *var, int dom,
                           const char * type, void *, DestructorFunction &df);

    void                   CalculateDomainNesting(void);

    void                   RegisterDataSelections(
                             const std::vector<avtDataSelection_p>&,
                             std::vector<bool>* applied
                           );
  protected:
    vtkAMAZEReader *myreader;
    bool           initializedReader;
    void           InitializeReader(void);
    bool           ReadShiftedGridInfo;
    bool           ComputeLog10;
    bool           ApplyLengthScale;
    ScaleOption    ScaleChoice;
    virtual void   PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual bool   ReturnsValidTime() const { return true; }
    virtual double GetTime(void);
    //virtual int    GetCycle(void);
    int            cycle;
    size_t         resolution; // for user selection of resolution
    virtual bool   HasInvariantMetaData(void) const { return false; };
    virtual bool   HasInvariantSIL(void) const { return false; };
    virtual int    GetCycleFromFilename(const char *f) const;
};


#endif
