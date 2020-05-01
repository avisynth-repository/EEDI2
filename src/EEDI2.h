/*
**                    EEDI2 v0.9.2 for AviSynth 2.5.x
**
**   EEDI2 resizes an image by 2x in the vertical direction by copying the
**   existing image to 2*y(n) and interpolating the missing field.  It is
**   intended for edge-directed interpolation for deinterlacing (i.e. not
**   really made for resizing a normal image).
**   
**   Copyright (C) 2005-2006 Kevin Stone
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <windows.h>
#include <math.h>
#include <malloc.h>
#include "Avisynth.h"
#include "PlanarFrame.h"

__declspec(align(16)) const int limlut[33] = { 
                         6, 6, 7, 7, 8, 8, 9, 9, 9, 10,
                         10, 11, 11, 12, 12, 12, 12, 12, 12, 12,
                         12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                         12, -1, -1 };

class EEDI2 : public GenericVideoFilter
{
private:
	int mthresh, vthresh, lthresh, pp, map;
	int estr, dstr, maxd, field, fieldS, nt;
	int *cx2, *cy2, *cxy, *tmpc;
	VideoInfo vi_saved;
	PlanarFrame srcPF, dstPF, mskPF, tmpPF, msk2PF;
	PlanarFrame dst2PF, dst2MPF, tmp2PF, tmp2PF2;
	void EEDI2::buildEdgeMask(int plane, PlanarFrame &msk, PlanarFrame &src);
	void EEDI2::dialate(int plane, PlanarFrame &msk, PlanarFrame &dst);
	void EEDI2::erode(int plane, PlanarFrame &msk, PlanarFrame &dst);
	void EEDI2::removeSmallHorzGaps(const int plane, PlanarFrame &msk, PlanarFrame &dst);
	void EEDI2::calcDirections(int plane, PlanarFrame &msk, PlanarFrame &src, PlanarFrame &dst);
	void EEDI2::filterDirMap(int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::filterDirMap2X(int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::expandDirMap(int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::expandDirMap2X(int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::upscaleBy2(int plane, PlanarFrame &src, PlanarFrame &dst, IScriptEnvironment *env);
	void EEDI2::markDirections2X(int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::InterpolateLattice(int plane, PlanarFrame &dmsk, PlanarFrame &dst,
		PlanarFrame &omsk, IScriptEnvironment *env);
	void EEDI2::sort_metrics(int *order, int length);
	void EEDI2::calcDerivatives(int plane, PlanarFrame &src, int *x2, int *y2, int *xy);
	void EEDI2::gaussianBlurSqrt2(int *src, int *tmp, int *dst, int width, int height, int pitch);
	void EEDI2::gaussianBlur1(int plane, PlanarFrame &src, PlanarFrame &tmp, PlanarFrame &dst);
	void EEDI2::postProcess(int plane, PlanarFrame &nmsk, PlanarFrame &omsk, PlanarFrame &src);
	void EEDI2::postProcessCorner(int plane, int *x2, int *y2, int *xy, PlanarFrame &msk, 
		PlanarFrame &dst, int pitch);
	void EEDI2::fillGaps2X(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);
	void EEDI2::filterMap(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst);

public:
	PVideoFrame __stdcall EEDI2::GetFrame(int n, IScriptEnvironment *env);
	EEDI2::EEDI2(PClip _child, int _mthresh, int _lthresh, int _vthresh, int _estr, int _dstr, 
		int _maxd, int _field, int _map, int _nt, int _pp, IScriptEnvironment *env);
	EEDI2::~EEDI2();
};