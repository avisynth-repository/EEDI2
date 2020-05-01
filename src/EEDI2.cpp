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

#include "EEDI2.h"

EEDI2::EEDI2(PClip _child, int _mthresh, int _lthresh, int _vthresh, int _estr, int _dstr, 
	int _maxd, int _field, int _map, int _nt, int _pp, IScriptEnvironment *env) : 
	GenericVideoFilter(_child), mthresh(_mthresh), lthresh(_lthresh), vthresh(_vthresh), 
	estr(_estr), dstr(_dstr), maxd(_maxd), field(_field), map(_map), nt(_nt), pp(_pp)
{
	cx2 = cy2 = cxy = tmpc = NULL;
	if (!vi.IsYUY2() && !vi.IsYV12())
		env->ThrowError("EEDI2:  YV12 or YUY2 input required!");
	if (field < -2 || field > 3)
		env->ThrowError("EEDI2:  field must be set to -2, -1, 0, 1, 2, or 3!");
	if (maxd >= 30)
		env->ThrowError("EEDI2:  maxd must be < 30!");
	if (map < 0 || map > 3)
		env->ThrowError("EEDI2:  map must be set to 0, 1, 2, or 3!");
	child->SetCacheHints(CACHE_NOTHING, 0);
	fieldS = field;
	if (field < 0) field = child->GetParity(0) ? 1 : 0;
	if (fieldS == -2) fieldS = field == 0 ? 2 : 3;
	else if (fieldS == -1) fieldS = field;
	else if (fieldS == 2) field = 0;
	else if (fieldS == 3) field = 1;
	srcPF.createFromProfile(vi);
	dstPF.createFromProfile(vi);
	mskPF.createFromProfile(vi);
	tmpPF.createFromProfile(vi);
	vi_saved = vi;
	if (map == 0 || map == 3)
	{
		vi.SetFieldBased(false);
		vi.height *= 2;
		dst2PF.createFromProfile(vi);
		dst2MPF.createFromProfile(vi);
		tmp2PF.createFromProfile(vi);
		tmp2PF2.createFromProfile(vi);
		msk2PF.createFromProfile(vi);
		if (pp > 1 && map == 0)
		{
			cx2 = (int*)_aligned_malloc(srcPF.GetHeight(0)*srcPF.GetPitch(0)*sizeof(int), 16);
			cy2 = (int*)_aligned_malloc(srcPF.GetHeight(0)*srcPF.GetPitch(0)*sizeof(int), 16);
			cxy = (int*)_aligned_malloc(srcPF.GetHeight(0)*srcPF.GetPitch(0)*sizeof(int), 16);
			tmpc = (int*)_aligned_malloc(srcPF.GetHeight(0)*srcPF.GetPitch(0)*sizeof(int), 16);
			if (!cx2 || !cy2 || !cxy || !tmpc)
				env->ThrowError("EEDI2:  malloc failure (pp>1)!");
		}
	}
	mthresh = mthresh*mthresh;
	vthresh = vthresh*81;
}

EEDI2::~EEDI2()
{
	if (cx2) _aligned_free(cx2);
	if (cy2) _aligned_free(cy2);
	if (cxy) _aligned_free(cxy);
	if (tmpc) _aligned_free(tmpc);
}

PVideoFrame __stdcall EEDI2::GetFrame(int n, IScriptEnvironment *env)
{
	if (fieldS > 1) field = (n&1) ? (fieldS == 2 ? 1 : 0) : (fieldS == 2 ? 0 : 1);
	srcPF.copyFrom(child->GetFrame(n,env),vi_saved);
	for (int b=0; b<3; ++b)
	{
		buildEdgeMask(b,mskPF,srcPF);
		erode(b,mskPF,tmpPF);
		dialate(b,tmpPF,mskPF);
		erode(b,mskPF,tmpPF);
		removeSmallHorzGaps(b,tmpPF,mskPF);
		if (map != 1)
		{
			calcDirections(b,mskPF,srcPF,tmpPF);
			filterDirMap(b,mskPF,tmpPF,dstPF);
			expandDirMap(b,mskPF,dstPF,tmpPF);
			filterMap(b,mskPF,tmpPF,dstPF);
			if (map == 2) continue;
			memset(dst2PF.GetPtr(b),0,dst2PF.GetHeight(b)*dst2PF.GetPitch(b));
			upscaleBy2(b,srcPF,dst2PF,env);
			upscaleBy2(b,dstPF,tmp2PF2,env);
			upscaleBy2(b,mskPF,msk2PF,env);
			markDirections2X(b,msk2PF,tmp2PF2,tmp2PF);
			filterDirMap2X(b,msk2PF,tmp2PF,dst2MPF);
			expandDirMap2X(b,msk2PF,dst2MPF,tmp2PF);
			fillGaps2X(b,msk2PF,tmp2PF,dst2MPF);
			fillGaps2X(b,msk2PF,dst2MPF,tmp2PF);
			if (map == 3) continue;
			InterpolateLattice(b,tmp2PF,dst2PF,tmp2PF2,env);
			if (pp)
			{
				if (pp == 1 || pp == 3)
				{
					tmp2PF.copyPlaneTo(tmp2PF2,b);
					filterDirMap2X(b,msk2PF,tmp2PF,dst2MPF);
					expandDirMap2X(b,msk2PF,dst2MPF,tmp2PF);
					postProcess(b,tmp2PF,tmp2PF2,dst2PF);
				}
				if (pp == 2 || pp == 3)
				{
					gaussianBlur1(b,srcPF,tmpPF,srcPF);
					calcDerivatives(b,srcPF,cx2,cy2,cxy);
					gaussianBlurSqrt2(cx2,tmpc,cx2,srcPF.GetWidth(b),srcPF.GetHeight(b),
						srcPF.GetPitch(b));
					gaussianBlurSqrt2(cy2,tmpc,cy2,srcPF.GetWidth(b),srcPF.GetHeight(b),
						srcPF.GetPitch(b));
					gaussianBlurSqrt2(cxy,tmpc,cxy,srcPF.GetWidth(b),srcPF.GetHeight(b),
						srcPF.GetPitch(b));
					postProcessCorner(b,cx2,cy2,cxy,tmp2PF2,dst2PF,srcPF.GetPitch(b));
				}
			}
		}
	}
	PVideoFrame dst = env->NewVideoFrame(vi);
	if (map == 0) dst2PF.copyTo(dst,vi);
	else if (map == 1) mskPF.copyTo(dst,vi);
	else if (map == 2) dstPF.copyTo(dst,vi);
	else tmp2PF.copyTo(dst,vi);
	return dst;
}

void EEDI2::buildEdgeMask(const int plane, PlanarFrame &msk, PlanarFrame &src)
{
	memset(msk.GetPtr(plane),0,msk.GetHeight(plane)*msk.GetPitch(plane));
	unsigned char *srcp = src.GetPtr(plane);
	unsigned char *dstp = msk.GetPtr(plane);
	const int height = src.GetHeight(plane);
	const int width = src.GetWidth(plane);
	const int src_pitch = src.GetPitch(plane);
	const int dst_pitch = msk.GetPitch(plane);
	srcp += src_pitch;
	dstp += dst_pitch;
	unsigned char *srcpp = srcp-src_pitch;
	unsigned char *srcpn = srcp+src_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if ((abs(srcpp[x]-srcp[x]) < 10 && abs(srcp[x]-srcpn[x]) < 10 &&
				 abs(srcpp[x]-srcpn[x]) < 10) ||
				(abs(srcpp[x-1]-srcp[x-1]) < 10 && abs(srcp[x-1]-srcpn[x-1]) < 10 &&
				 abs(srcpp[x-1]-srcpn[x-1]) < 10 && abs(srcpp[x+1]-srcp[x+1]) < 10 && 
				 abs(srcp[x+1]-srcpn[x+1]) < 10 && abs(srcpp[x+1]-srcpn[x+1]) < 10))
					continue;
			const int sum = srcpp[x-1]+srcpp[x]+srcpp[x+1]+
				srcp[x-1]+srcp[x]+srcp[x+1]+
				srcpn[x-1]+srcpn[x]+srcpn[x+1];
			const int sumsq = srcpp[x-1]*srcpp[x-1]+srcpp[x]*srcpp[x]+srcpp[x+1]*srcpp[x+1]+
				srcp[x-1]*srcp[x-1]+srcp[x]*srcp[x]+srcp[x+1]*srcp[x+1]+
				srcpn[x-1]*srcpn[x-1]+srcpn[x]*srcpn[x]+srcpn[x+1]*srcpn[x+1];
			if (9*sumsq-sum*sum < vthresh) continue;
			const int Ix = srcp[x+1]-srcp[x-1];
			const int Iy = max(max(abs(srcpp[x]-srcpn[x]),abs(srcpp[x]-srcp[x])),abs(srcp[x]-srcpn[x]));
			if (Ix*Ix+Iy*Iy >= mthresh)
			{
				dstp[x] = 255;
				continue;
			}
			const int Ixx = srcp[x-1]-2*srcp[x]+srcp[x+1];
			const int Iyy = srcpp[x]-2*srcp[x]+srcpn[x];
			if (abs(Ixx)+abs(Iyy) >= lthresh) dstp[x] = 255;
		}
		dstp += dst_pitch;
		srcpp += src_pitch;
		srcp += src_pitch;
		srcpn += src_pitch;
	}
}

void EEDI2::dialate(const int plane, PlanarFrame &msk, PlanarFrame &dst)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,mskp,msk_pitch,width,height);
	mskp += msk_pitch;
	unsigned char *mskpp = mskp-msk_pitch;
	unsigned char *mskpn = mskp+msk_pitch;
	dstp += dst_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0) continue;
			int count = 0;
			if (mskpp[x-1] == 0xFF) ++count;
			if (mskpp[x] == 0xFF) ++count;
			if (mskpp[x+1] == 0xFF) ++count;
			if (mskp[x-1] == 0xFF) ++count;
			if (mskp[x+1] == 0xFF) ++count;
			if (mskpn[x-1] == 0xFF) ++count;
			if (mskpn[x] == 0xFF) ++count;
			if (mskpn[x+1] == 0xFF) ++count;
			if (count >= dstr) dstp[x] = 0xFF;
		}
		mskpp += msk_pitch;
		mskp += msk_pitch;
		mskpn += msk_pitch;
		dstp += dst_pitch;
	}
}

void EEDI2::erode(const int plane, PlanarFrame &msk, PlanarFrame &dst)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,mskp,msk_pitch,width,height);
	mskp += msk_pitch;
	unsigned char *mskpp = mskp-msk_pitch;
	unsigned char *mskpn = mskp+msk_pitch;
	dstp += dst_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0xFF) continue;
			int count = 0;
			if (mskpp[x-1] == 0xFF) ++count;
			if (mskpp[x] == 0xFF) ++count;
			if (mskpp[x+1] == 0xFF) ++count;
			if (mskp[x-1] == 0xFF) ++count;
			if (mskp[x+1] == 0xFF) ++count;
			if (mskpn[x-1] == 0xFF) ++count;
			if (mskpn[x] == 0xFF) ++count;
			if (mskpn[x+1] == 0xFF) ++count;
			if (count < estr) dstp[x] = 0;
		}
		mskpp += msk_pitch;
		mskp += msk_pitch;
		mskpn += msk_pitch;
		dstp += dst_pitch;
	}
}

void EEDI2::removeSmallHorzGaps(const int plane, PlanarFrame &msk, PlanarFrame &dst)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,mskp,msk_pitch,width,height);
	mskp += msk_pitch;
	dstp += dst_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=3; x<width-3; ++x)
		{
			if (mskp[x])
			{
				if (mskp[x-3]) continue;
				if (mskp[x-2]) continue;
				if (mskp[x-1]) continue;
				if (mskp[x+1]) continue;
				if (mskp[x+2]) continue;
				if (mskp[x+3]) continue;
				dstp[x] = 0;
			}
			else
			{
				if ((mskp[x+1] && (mskp[x-1] || mskp[x-2] || mskp[x-3])) ||
					(mskp[x+2] && (mskp[x-1] || mskp[x-2])) ||
					(mskp[x+3] && mskp[x-1]))
					dstp[x] = 0xFF;
			}
		}
		mskp += msk_pitch;
		dstp += dst_pitch;
	}
}

void EEDI2::calcDirections(const int plane, PlanarFrame &msk, PlanarFrame &src, PlanarFrame &dst)
{
	memset(dst.GetPtr(plane),255,dst.GetPitch(plane)*dst.GetHeight(plane));
	unsigned char *mskp = msk.GetPtr(plane);
	unsigned char *srcp = src.GetPtr(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int src_pitch = src.GetPitch(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int width = src.GetWidth(plane);
	const int height = src.GetHeight(plane);
	mskp += msk_pitch;
	dstp += dst_pitch;
	srcp += src_pitch;
	unsigned char *src2p = srcp-src_pitch*2;
	unsigned char *srcpp = srcp-src_pitch;
	unsigned char *srcpn = srcp+src_pitch;
	unsigned char *src2n = srcp+src_pitch*2;
	unsigned char *mskpp = mskp-msk_pitch;
	unsigned char *mskpn = mskp+msk_pitch;
	const int maxdt = plane == 0 ? maxd : (maxd>>1);
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0xFF || (mskp[x-1] != 0xFF && mskp[x+1] != 0xFF)) 
				continue;
			const int startu = max(-x+1,-maxdt);
			const int stopu = min(width-2-x,maxdt);
			int minb = min(13*nt,(abs(srcp[x]-srcpn[x])+abs(srcp[x]-srcpp[x]))*6);
			int mina = min(19*nt,(abs(srcp[x]-srcpn[x])+abs(srcp[x]-srcpp[x]))*9);
			int minc = mina;
			int mind = minb;
			int mine = minb;
			int dira = -5000, dirb = -5000, dirc = -5000, dird = -5000, dire = -5000;
			for (int u=startu; u<=stopu; ++u)
			{
				if ((y == 1 || mskpp[x-1+u] == 0xFF || mskpp[x+u] == 0xFF || mskpp[x+1+u] == 0xFF) &&
					(y == height-2 || mskpn[x-1-u] == 0xFF || mskpn[x-u] == 0xFF || mskpn[x+1-u] == 0xFF))
				{
					const int diffsn = abs(srcp[x-1]-srcpn[x-1-u])+abs(srcp[x]-srcpn[x-u])+abs(srcp[x+1]-srcpn[x+1-u]);
					const int diffsp = abs(srcp[x-1]-srcpp[x-1+u])+abs(srcp[x]-srcpp[x+u])+abs(srcp[x+1]-srcpp[x+1+u]);
					const int diffps = abs(srcpp[x-1]-srcp[x-1-u])+abs(srcpp[x]-srcp[x-u])+abs(srcpp[x+1]-srcp[x+1-u]);
					const int diffns = abs(srcpn[x-1]-srcp[x-1+u])+abs(srcpn[x]-srcp[x+u])+abs(srcpn[x+1]-srcp[x+1+u]);
					const int diff = diffsn+diffsp+diffps+diffns;
					int diffd = diffsp+diffns;
					int diffe = diffsn+diffps;
					if (diff < minb)
					{
						dirb = u;
						minb = diff;
					}
					if (y > 1)
					{
						const diff2pp = abs(src2p[x-1]-srcpp[x-1-u])+abs(src2p[x]-srcpp[x-u])+abs(src2p[x+1]-srcpp[x+1-u]);
						const diffp2p = abs(srcpp[x-1]-src2p[x-1+u])+abs(srcpp[x]-src2p[x+u])+abs(srcpp[x+1]-src2p[x+1+u]);
						const int diffa = diff+diff2pp+diffp2p;
						diffd += diffp2p;
						diffe += diff2pp;
						if (diffa < mina)
						{
							dira = u;
							mina = diffa;
						}
					}
					if (y < height-2)
					{
						const int diff2nn = abs(src2n[x-1]-srcpn[x-1+u])+abs(src2n[x]-srcpn[x+u])+abs(src2n[x+1]-srcpn[x+1+u]);
						const int diffn2n = abs(srcpn[x-1]-src2n[x-1-u])+abs(srcpn[x]-src2n[x-u])+abs(srcpn[x+1]-src2n[x+1-u]);
						const int diffc = diff+diff2nn+diffn2n;
						diffd += diff2nn;
						diffe += diffn2n;
						if (diffc < minc)
						{
							dirc = u;
							minc = diffc;
						}
					}
					if (diffd < mind)
					{
						dird = u;
						mind = diffd;
					}
					if (diffe < mine)
					{
						dire = u;
						mine = diffe;
					}
				}
			}
			int order[5], k=0;
			if (dira != -5000) order[k++] = dira;
			if (dirb != -5000) order[k++] = dirb;
			if (dirc != -5000) order[k++] = dirc;
			if (dird != -5000) order[k++] = dird;
			if (dire != -5000) order[k++] = dire;
			if (k > 1)
			{
				sort_metrics(order, k);
				const int mid = (k&1) ? order[k>>1] : (order[(k-1)>>1]+order[k>>1]+1)>>1;
				const int tlim = max(limlut[abs(mid)]>>2,2);
				int sum = 0, count = 0;
				for (int i=0; i<k; ++i)
				{
					if (abs(order[i]-mid) <= tlim)
					{
						++count;
						sum += order[i];
					}
				}
				if (count > 1) 
					dstp[x] = 128+(int(float(sum)/float(count)))*4;
				else dstp[x] = 128;
			}
			else dstp[x] = 128;
		}
		mskpp += msk_pitch;
		mskp += msk_pitch;
		mskpn += msk_pitch;
		src2p += src_pitch;
		srcpp += src_pitch;
		srcp += src_pitch;
		srcpn += src_pitch;
		src2n += src_pitch;
		dstp += dst_pitch;
	}
}

void EEDI2::filterMap(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	mskp += msk_pitch;
	dmskp += dmsk_pitch;
	dstp += dst_pitch;
	unsigned char *dmskpp = dmskp-dmsk_pitch;
	unsigned char *dmskpn = dmskp+dmsk_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (dmskp[x] == 0xFF || mskp[x] != 0xFF) continue;
			const int dir = (dmskp[x]-128)>>2;
			const int lim = max(abs(dir)*2,12);
			bool ict = false, icb = false;
			if (dir < 0)
			{
				const int dirt = max(-x,dir);
				for (int j=dirt; j<=0; ++j)
				{
					if ((abs(dmskpp[x+j]-dmskp[x]) > lim && dmskpp[x+j] != 0xFF) ||
						(dmskp[x+j] == 0xFF && dmskpp[x+j] == 0xFF) ||
						(abs(dmskp[x+j]-dmskp[x]) > lim && dmskp[x+j] != 0xFF))
					{
						ict = true;
						break;
					}
				}
			}
			else
			{
				const int dirt = min(width-x-1,dir);
				for (int j=0; j<=dirt; ++j)
				{
					if ((abs(dmskpp[x+j]-dmskp[x]) > lim && dmskpp[x+j] != 0xFF) ||
						(dmskp[x+j] == 0xFF && dmskpp[x+j] == 0xFF) ||
						(abs(dmskp[x+j]-dmskp[x]) > lim && dmskp[x+j] != 0xFF))
					{
						ict = true;
						break;
					}
				}
			}
			if (ict)
			{
				if (dir < 0)
				{
					const int dirt = min(width-x-1,abs(dir));
					for (int j=0; j<=dirt; ++j)
					{
						if ((abs(dmskpn[x+j]-dmskp[x]) > lim && dmskpn[x+j] != 0xFF) ||
							(dmskpn[x+j] == 0xFF && dmskp[x+j] == 0xFF) ||
							(abs(dmskp[x+j]-dmskp[x]) > lim && dmskp[x+j] != 0xFF))
						{
							icb = true;
							break;
						}
					}
				}
				else
				{
					const int dirt = max(-x,-dir);
					for (int j=dirt; j<=0; ++j)
					{
						if ((abs(dmskpn[x+j]-dmskp[x]) > lim && dmskpn[x+j] != 0xFF) ||
							(dmskpn[x+j] == 0xFF && dmskp[x+j] == 0xFF) ||
							(abs(dmskp[x+j]-dmskp[x]) > lim && dmskp[x+j] != 0xFF))
						{
							icb = true;
							break;
						}
					}
				}
				if (icb) dstp[x] = 255;
			}
		}
		mskp += msk_pitch;
		dmskpp += dmsk_pitch;
		dmskp += dmsk_pitch;
		dmskpn += dmsk_pitch;
		dstp += dst_pitch;
	}
}

void EEDI2::filterDirMap(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	dmskp += dmsk_pitch;
	unsigned char *dmskpp = dmskp-dmsk_pitch;
	unsigned char *dmskpn = dmskp+dmsk_pitch;
	dstp += dst_pitch;
	mskp += msk_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0xFF) continue;
			int u = 0, order[9];
			if (dmskpp[x-1] != 0xFF) order[u++] = dmskpp[x-1];
			if (dmskpp[x] != 0xFF) order[u++] = dmskpp[x];
			if (dmskpp[x+1] != 0xFF) order[u++] = dmskpp[x+1];
			if (dmskp[x-1] != 0xFF) order[u++] = dmskp[x-1];
			if (dmskp[x] != 0xFF) order[u++] = dmskp[x];
			if (dmskp[x+1] != 0xFF) order[u++] = dmskp[x+1];
			if (dmskpn[x-1] != 0xFF) order[u++] = dmskpn[x-1];
			if (dmskpn[x] != 0xFF) order[u++] = dmskpn[x];
			if (dmskpn[x+1] != 0xFF) order[u++] = dmskpn[x+1];
			if (u < 4)
			{
				dstp[x] = 255;
				continue;
			}
			sort_metrics(order, u);
			const int mid = (u&1) ? order[u>>1] : (order[(u-1)>>1]+order[u>>1]+1)>>1;
			int sum = 0, count = 0;
			const int lim = limlut[abs(mid-128)>>2];
			for (int i=0; i<u; ++i)
			{
				if (abs(order[i]-mid) <= lim)
				{
					++count;
					sum += order[i];
				}
			}
			if (count < 4 || (count < 5 && dmskp[x] == 0xFF))
			{
				dstp[x] = 255;
				continue;
			}
			dstp[x] = int((float(sum+mid)/float(count+1))+0.5f);
		}
		dmskpp += dmsk_pitch;
		dmskp += dmsk_pitch;
		dmskpn += dmsk_pitch;
		dstp += dst_pitch;
		mskp += msk_pitch;
	}
}

void EEDI2::filterDirMap2X(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	dmskp += dmsk_pitch*(2-field);
	unsigned char *dmskpp = dmskp-dmsk_pitch*2;
	unsigned char *dmskpn = dmskp+dmsk_pitch*2;
	mskp += msk_pitch*(1-field);
	unsigned char *mskpn = mskp+msk_pitch*2;
	dstp += dst_pitch*(2-field);
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0xFF && mskpn[x] != 0xFF) continue;
			int u = 0, order[9];
			if (y > 1)
			{
				if (dmskpp[x-1] != 0xFF) order[u++] = dmskpp[x-1];
				if (dmskpp[x] != 0xFF) order[u++] = dmskpp[x];
				if (dmskpp[x+1] != 0xFF) order[u++] = dmskpp[x+1];
			}
			if (dmskp[x-1] != 0xFF) order[u++] = dmskp[x-1];
			if (dmskp[x] != 0xFF) order[u++] = dmskp[x];
			if (dmskp[x+1] != 0xFF) order[u++] = dmskp[x+1];
			if (y < height-2)
			{
				if (dmskpn[x-1] != 0xFF) order[u++] = dmskpn[x-1];
				if (dmskpn[x] != 0xFF) order[u++] = dmskpn[x];
				if (dmskpn[x+1] != 0xFF) order[u++] = dmskpn[x+1];
			}
			if (u < 4)
			{
				dstp[x] = 255;
				continue;
			}
			sort_metrics(order, u);
			const int mid = (u&1) ? order[u>>1] : (order[(u-1)>>1]+order[u>>1]+1)>>1;
			int sum = 0, count = 0;
			const int lim = limlut[abs(mid-128)>>2];
			for (int i=0; i<u; ++i)
			{
				if (abs(order[i]-mid) <= lim)
				{
					++count;
					sum += order[i];
				}
			}
			if (count < 4 || (count < 5 && dmskp[x] == 0xFF))
			{
				dstp[x] = 255;
				continue;
			}
			dstp[x] = int((float(sum+mid)/float(count+1))+0.5f);
		}
		mskp += msk_pitch*2;
		mskpn += msk_pitch*2;
		dmskpp += dmsk_pitch*2;
		dmskp += dmsk_pitch*2;
		dmskpn += dmsk_pitch*2;
		dstp += dst_pitch*2;
	}
}

void EEDI2::expandDirMap(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	dmskp += dmsk_pitch;
	unsigned char *dmskpp = dmskp-dmsk_pitch;
	unsigned char *dmskpn = dmskp+dmsk_pitch;
	dstp += dst_pitch;
	mskp += msk_pitch;
	for (int y=1; y<height-1; ++y)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (dmskp[x] != 0xFF || mskp[x] != 0xFF) continue;
			int u = 0, order[9];
			if (dmskpp[x-1] != 0xFF) order[u++] = dmskpp[x-1];
			if (dmskpp[x] != 0xFF) order[u++] = dmskpp[x];
			if (dmskpp[x+1] != 0xFF) order[u++] = dmskpp[x+1];
			if (dmskp[x-1] != 0xFF) order[u++] = dmskp[x-1];
			if (dmskp[x+1] != 0xFF) order[u++] = dmskp[x+1];
			if (dmskpn[x-1] != 0xFF) order[u++] = dmskpn[x-1];
			if (dmskpn[x] != 0xFF) order[u++] = dmskpn[x];
			if (dmskpn[x+1] != 0xFF) order[u++] = dmskpn[x+1];
			if (u < 5) continue;
			sort_metrics(order, u);
			const int mid = (u&1) ? order[u>>1] : (order[(u-1)>>1]+order[u>>1]+1)>>1;
			int sum = 0, count = 0;
			const int lim = limlut[abs(mid-128)>>2];
			for (int i=0; i<u; ++i)
			{
				if (abs(order[i]-mid) <= lim)
				{
					++count;
					sum += order[i];
				}
			}
			if (count < 5) continue;
			dstp[x] = int((float(sum+mid)/float(count+1))+0.5f);
		}
		dmskpp += dmsk_pitch;
		dmskp += dmsk_pitch;
		dmskpn += dmsk_pitch;
		dstp += dst_pitch;
		mskp += msk_pitch;
	}
}

void EEDI2::expandDirMap2X(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	dmskp += dmsk_pitch*(2-field);
	unsigned char *dmskpp = dmskp-dmsk_pitch*2;
	unsigned char *dmskpn = dmskp+dmsk_pitch*2;
	mskp += msk_pitch*(1-field);
	unsigned char *mskpn = mskp+msk_pitch*2;
	dstp += dst_pitch*(2-field);
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (dmskp[x] != 0xFF || (mskp[x] != 0xFF && mskpn[x] != 0xFF)) continue;
			int u = 0, order[9];
			if (y > 1)
			{
				if (dmskpp[x-1] != 0xFF) order[u++] = dmskpp[x-1];
				if (dmskpp[x] != 0xFF) order[u++] = dmskpp[x];
				if (dmskpp[x+1] != 0xFF) order[u++] = dmskpp[x+1];
			}
			if (dmskp[x-1] != 0xFF) order[u++] = dmskp[x-1];
			if (dmskp[x+1] != 0xFF) order[u++] = dmskp[x+1];
			if (y < height-2)
			{
				if (dmskpn[x-1] != 0xFF) order[u++] = dmskpn[x-1];
				if (dmskpn[x] != 0xFF) order[u++] = dmskpn[x];
				if (dmskpn[x+1] != 0xFF) order[u++] = dmskpn[x+1];
			}
			if (u < 5) continue;
			sort_metrics(order, u);
			const int mid = (u&1) ? order[u>>1] : (order[(u-1)>>1]+order[u>>1]+1)>>1;
			int sum = 0, count = 0;
			const int lim = limlut[abs(mid-128)>>2];
			for (int i=0; i<u; ++i)
			{
				if (abs(order[i]-mid) <= lim)
				{
					++count;
					sum += order[i];
				}
			}
			if (count < 5) continue;
			dstp[x] = int((float(sum+mid)/float(count+1))+0.5f);
		}
		mskp += msk_pitch*2;
		mskpn += msk_pitch*2;
		dmskpp += dmsk_pitch*2;
		dmskp += dmsk_pitch*2;
		dmskpn += dmsk_pitch*2;
		dstp += dst_pitch*2;
	}
}

void EEDI2::fillGaps2X(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	unsigned char *dmskp = dmsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	dst.BitBlt(dstp,dst_pitch,dmskp,dmsk_pitch,width,height);
	dmskp += dmsk_pitch*(2-field);
	unsigned char *dmskpp = dmskp-dmsk_pitch*2;
	unsigned char *dmskpn = dmskp+dmsk_pitch*2;
	mskp += msk_pitch*(1-field);
	unsigned char *mskpp = mskp-msk_pitch*2;
	unsigned char *mskpn = mskp+msk_pitch*2;
	unsigned char *mskpnn = mskpn+msk_pitch*2;
	dstp += dst_pitch*(2-field);
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (dmskp[x] != 0xFF || 
				(mskp[x] != 0xFF && mskpn[x] != 0xFF)) continue;
			int u = x-1, back = 500, forward = -500;
			while(u)
			{
				if (dmskp[u] != 0xFF) 
				{ 
					back = dmskp[u]; 
					break; 
				}
				if (mskp[u] != 0xFF && mskpn[u] != 0xFF) break;
				--u;
			}
			int v = x+1;
			while (v < width)
			{
				if (dmskp[v] != 0xFF)
				{
					forward = dmskp[v];
					break;
				}
				if (mskp[v] != 0xFF && mskpn[v] != 0xFF) break;
				++v;
			}
			bool tc = true, bc = true;
			int mint = 500, maxt = -20;
			int minb = 500, maxb = -20;
			for (int j=u; j<=v; ++j)
			{
				if (tc)
				{
					if (y <= 2 || dmskpp[j] == 0xFF || (mskpp[j] != 0xFF && mskp[j] != 0xFF))
					{
						tc = false;
						mint = maxt = 20;
					}
					else
					{
						if (dmskpp[j] < mint) mint = dmskpp[j];
						if (dmskpp[j] > maxt) maxt = dmskpp[j];
					}
				}
				if (bc)
				{
					if (y >= height-3 || dmskpn[j] == 0xFF || (mskpn[j] != 0xFF && mskpnn[j] != 0xFF))
					{
						bc = false;
						minb = maxb = 20;
					}
					else
					{
						if (dmskpn[j] < minb) minb = dmskpn[j];
						if (dmskpn[j] > maxb) maxb = dmskpn[j];
					}
				}
			}
			if (maxt == -20) maxt = mint = 20;
			if (maxb == -20) maxb = minb = 20;
			int thresh = max(max(max(abs(forward-128),abs(back-128))>>2,8),max(abs(mint-maxt),abs(minb-maxb)));
			const int flim = min(max(abs(forward-128),abs(back-128))>>2,6);
			if (abs(forward-back) <= thresh && (v-u-1 <= flim || tc || bc))
			{
				double step = double(forward-back)/double(v-u);
				for (int j=0; j<v-u-1; ++j)
					dstp[u+j+1] = back+int(j*step+0.5);
			}
		}
		mskpp += msk_pitch*2;
		mskp += msk_pitch*2;
		mskpn += msk_pitch*2;
		mskpnn += msk_pitch*2;
		dmskpp += dmsk_pitch*2;
		dmskp += dmsk_pitch*2;
		dmskpn += dmsk_pitch*2;
		dstp += dst_pitch*2;
	}
}

void EEDI2::sort_metrics(int *order, const int length)
{
	for (int i=1; i<length; ++i) 
	{
		int j = i;
		const int temp = order[j];
		while (j>0 && order[j-1]>temp) 
		{
			order[j] = order[j-1];
			--j;
		}
		order[j] = temp;
	}
}

void EEDI2::upscaleBy2(const int plane, PlanarFrame &src, PlanarFrame &dst, IScriptEnvironment *env)
{
	env->BitBlt(dst.GetPtr(plane)+(1-field)*dst.GetPitch(plane),dst.GetPitch(plane)*2,src.GetPtr(plane),
		src.GetPitch(plane),src.GetWidth(plane),src.GetHeight(plane));
}

void EEDI2::markDirections2X(const int plane, PlanarFrame &msk, PlanarFrame &dmsk, PlanarFrame &dst)
{
	unsigned char *dmskp = dmsk.GetPtr(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	unsigned char *mskp = msk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	memset(dstp,255,dst_pitch*height);
	dstp += dst_pitch*(2-field);
	dmskp += dmsk_pitch*(1-field);
	mskp += msk_pitch*(1-field);
	unsigned char *dmskpn = dmskp+dmsk_pitch*2;
	unsigned char *mskpn = mskp+msk_pitch*2;
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=1; x<width-1; ++x)
		{
			if (mskp[x] != 0xFF && mskpn[x] != 0xFF) continue;
			int v = 0, order[6];
			if (dmskp[x-1] != 0xFF) order[v++] = dmskp[x-1];
			if (dmskp[x] != 0xFF) order[v++] = dmskp[x];
			if (dmskp[x+1] != 0xFF) order[v++] = dmskp[x+1];
			if (dmskpn[x-1] != 0xFF) order[v++] = dmskpn[x-1];
			if (dmskpn[x] != 0xFF) order[v++] = dmskpn[x];
			if (dmskpn[x+1] != 0xFF) order[v++] = dmskpn[x+1];
			if (v < 3) continue;
			else
			{
				sort_metrics(order,v);
				const int mid = (v&1) ? order[v>>1] : (order[(v-1)>>1]+order[v>>1]+1)>>1;
				const int lim = limlut[abs(mid-128)>>2];
				int u = 0;
				if (abs(dmskp[x-1]-dmskpn[x-1]) <= lim || dmskp[x-1] == 0xFF || dmskpn[x-1] == 0xFF) ++u;
				if (abs(dmskp[x]-dmskpn[x]) <= lim || dmskp[x] == 0xFF || dmskpn[x] == 0xFF) ++u;
				if (abs(dmskp[x+1]-dmskpn[x-1]) <= lim || dmskp[x+1] == 0xFF || dmskpn[x+1] == 0xFF) ++u;
				if (u < 2) continue;
				int count = 0, sum = 0;
				for (int i=0; i<v; ++i)
				{
					if (abs(order[i]-mid) <= lim)
					{
						++count;
						sum += order[i];
					}
				}
				if (count < v-2 || count < 2) continue;
				dstp[x] = int((float(sum+mid)/float(count+1))+0.5f);
			}
		}
		mskp += msk_pitch*2;
		mskpn += msk_pitch*2;
		dstp += dst_pitch*2;
		dmskp += dmsk_pitch*2;
		dmskpn += dmsk_pitch*2;
	}
}

void EEDI2::InterpolateLattice(const int plane, PlanarFrame &dmsk, PlanarFrame &dst,
							  PlanarFrame &omsk, IScriptEnvironment *env)
{
	unsigned char *dmskp = dmsk.GetPtr(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	unsigned char *omskp = omsk.GetPtr(plane);
	const int dmsk_pitch = dmsk.GetPitch(plane);
	const int dst_pitch = dst.GetPitch(plane);
	const int omsk_pitch = omsk.GetPitch(plane);
	const int height = dst.GetHeight(plane);
	const int width = dst.GetWidth(plane);
	if (field == 1) env->BitBlt(dstp+(height-1)*dst_pitch,dst_pitch,dstp+(height-2)*dst_pitch,dst_pitch,width,1);
	else env->BitBlt(dstp,dst_pitch,dstp+dst_pitch,dst_pitch,width,1);
	dstp += dst_pitch*(1-field);
	omskp += omsk_pitch*(1-field);
	unsigned char *dstpn = dstp+dst_pitch;
	unsigned char *dstpnn = dstp+dst_pitch*2;
	unsigned char *omskn = omskp + omsk_pitch*2;
	dmskp += dmsk_pitch*(2-field);
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=0; x<width; ++x)
		{
			int dir = dmskp[x];
			const int lim = limlut[abs(dir-128)>>2];
			if (dir == 255 || (abs(dmskp[x]-dmskp[x-1]) > lim && abs(dmskp[x]-dmskp[x+1]) > lim))
			{
				dstpn[x] = (dstp[x]+dstpnn[x]+1)>>1;
				if (dir != 255) dmskp[x] = 128;
				continue;
			}
			if (lim < 9)
			{
				const int sum = dstp[x-1]+dstp[x]+dstp[x+1]+
					dstpnn[x-1]+dstpnn[x]+dstpnn[x+1];
				const int sumsq = dstp[x-1]*dstp[x-1]+dstp[x]*dstp[x]+dstp[x+1]*dstp[x+1]+
					dstpnn[x-1]*dstpnn[x-1]+dstpnn[x]*dstpnn[x]+dstpnn[x+1]*dstpnn[x+1];
				if (6*sumsq-sum*sum < 576)
				{
					dstpn[x] = (dstp[x]+dstpnn[x]+1)>>1;
					dmskp[x] = 255;
					continue;
				}
			}
			if (x > 1 && x < width-2 && 
				(dstp[x] < max(dstp[x-2],dstp[x-1])-3 && dstp[x] < max(dstp[x+2],dstp[x+1])-3 &&
				dstpnn[x] < max(dstpnn[x-2],dstpnn[x-1])-3 && dstpnn[x] < max(dstpnn[x+2],dstpnn[x+1])-3) ||
				(dstp[x] > min(dstp[x-2],dstp[x-1])+3 && dstp[x] > min(dstp[x+2],dstp[x+1])+3 &&
				dstpnn[x]  > min(dstpnn[x-2],dstpnn[x-1])+3 && dstpnn[x] > min(dstpnn[x+2],dstpnn[x+1])+3))
			{
				dstpn[x] = (dstp[x]+dstpnn[x]+1)>>1;
				dmskp[x] = 128;
				continue;
			}
			dir = (dir-128+2)>>2;
			int val = (dstp[x]+dstpnn[x]+1)>>1;
			const int startu = (dir-2 < 0) ? max(-x+1,max(dir-2,-width+2+x)) :
				min(x-1,min(dir-2,width-2-x));
			const int stopu = (dir+2 < 0) ? max(-x+1,max(dir+2,-width+2+x)) :
				min(x-1,min(dir+2,width-2-x));
			int min = 8*nt;
			for (int u=startu; u<=stopu; ++u)
			{
				const int diff = abs(dstp[x-1]-dstpnn[x-u-1])+abs(dstp[x]-dstpnn[x-u])+abs(dstp[x+1]-dstpnn[x-u+1])+
						abs(dstpnn[x-1]-dstp[x+u-1])+abs(dstpnn[x]-dstp[x+u])+abs(dstpnn[x+1]-dstp[x+u+1]);
				if (diff < min && 
					((omskp[x-1+u] != 0xFF && abs(omskp[x-1+u]-dmskp[x]) <= lim) ||
					 (omskp[x+u] != 0xFF && abs(omskp[x+u]-dmskp[x]) <= lim) ||
					 (omskp[x+1+u] != 0xFF && abs(omskp[x+1+u]-dmskp[x]) <= lim)) &&
					((omskn[x-1-u] != 0xFF && abs(omskn[x-1-u]-dmskp[x]) <= lim) ||
					 (omskn[x-u] != 0xFF && abs(omskn[x-u]-dmskp[x]) <= lim) ||
					 (omskn[x+1-u] != 0xFF && abs(omskn[x+1-u]-dmskp[x]) <= lim)))
				{
					const int diff2 = abs(dstp[x+(u>>1)-1]-dstpnn[x-(u>>1)-1])+abs(dstp[x+(u>>1)]-dstpnn[x-(u>>1)])+
							abs(dstp[x+(u>>1)+1]-dstpnn[x-(u>>1)+1]);
					if (diff2 < 4*nt && (((abs(omskp[x+(u>>1)]-omskn[x-(u>>1)]) <= lim ||
										   abs(omskp[x+(u>>1)]-omskn[x-((u+1)>>1)]) <= lim) && 
										   omskp[x+(u>>1)] != 0xFF) || 
										 ((abs(omskp[x+((u+1)>>1)]-omskn[x-(u>>1)]) <= lim ||
										   abs(omskp[x+((u+1)>>1)]-omskn[x-((u+1)>>1)]) <= lim) && 
										   omskp[x+((u+1)>>1)] != 0xFF))) 
					{
						if ((abs(dmskp[x]-omskp[x+(u>>1)]) <= lim || abs(dmskp[x]-omskp[x+((u+1)>>1)]) <= lim) &&
							(abs(dmskp[x]-omskn[x-(u>>1)]) <= lim || abs(dmskp[x]-omskn[x-((u+1)>>1)]) <= lim))
						{
							val = (dstp[x+(u>>1)]+dstp[x+((u+1)>>1)]+
								dstpnn[x-(u>>1)]+dstpnn[x-((u+1)>>1)]+2)>>2;
							min = diff;
							dir = u;
						}
					}
				}
			}
			if (min != 8*nt)
			{
				dstpn[x] = val;
				dmskp[x] = 128+dir*4;
			}
			else 
			{
				const int minm = min(dstp[x],dstpnn[x]);
				const int maxm = max(dstp[x],dstpnn[x]);
				const int d = plane == 0 ? 4 : 2;
				const int startu = max(-x+1,-d);
				const int stopu = min(width-2-x,d);
				min = 7*nt;
				for (int u=startu; u<=stopu; ++u)
				{
					const int p1 = dstp[x+(u>>1)]+dstp[x+((u+1)>>1)];
					const int p2 = dstpnn[x-(u>>1)]+dstpnn[x-((u+1)>>1)];
					const int diff = abs(dstp[x-1]-dstpnn[x-u-1])+abs(dstp[x]-dstpnn[x-u])+abs(dstp[x+1]-dstpnn[x-u+1])+
							abs(dstpnn[x-1]-dstp[x+u-1])+abs(dstpnn[x]-dstp[x+u])+abs(dstpnn[x+1]-dstp[x+u+1])+
							abs(p1-p2);
					if (diff < min)
					{
						const int valt = (p1+p2+2)>>2;
						if (valt >= minm && valt <= maxm)
						{
							val = valt;
							min = diff;
							dir = u;
						}
					}
				}
				dstpn[x] = val;
				if (min == 7*nt) dmskp[x] = 128;
				else dmskp[x] = 128+dir*4;
			}
		}
		dstp += dst_pitch*2;
		dstpn += dst_pitch*2;
		dstpnn += dst_pitch*2;
		dmskp += dmsk_pitch*2;
		omskp += omsk_pitch*2;
		omskn += omsk_pitch*2;
	}
}

void EEDI2::postProcess(const int plane, PlanarFrame &nmsk, PlanarFrame &omsk, PlanarFrame &src)
{
	unsigned char *nmskp = nmsk.GetPtr(plane);
	const int width = nmsk.GetWidth(plane);
	const int height = nmsk.GetHeight(plane);
	const int nmsk_pitch = nmsk.GetPitch(plane);
	unsigned char *omskp = omsk.GetPtr(plane);
	const int omsk_pitch = omsk.GetPitch(plane);
	unsigned char *dstp = src.GetPtr(plane);
	const int src_pitch = src.GetPitch(plane);
	nmskp += (2-field)*nmsk_pitch;
	omskp += (2-field)*omsk_pitch;
	dstp += (2-field)*src_pitch;
	unsigned char *srcpp = dstp-src_pitch;
	unsigned char *srcpn = dstp+src_pitch;
	for (int y=2-field; y<height-1; y+=2)
	{
		for (int x=0; x<width; ++x)
		{
			const int lim = limlut[abs(nmskp[x]-128)>>2];
			if (abs(nmskp[x]-omskp[x]) > lim && omskp[x] != 255 && omskp[x] != 128)
				dstp[x] = (srcpp[x]+srcpn[x]+1)>>1;
		}
		nmskp += nmsk_pitch*2;
		omskp += omsk_pitch*2;
		srcpp += src_pitch*2;
		dstp += src_pitch*2;
		srcpn += src_pitch*2;
	}
}

void EEDI2::postProcessCorner(const int plane, int *x2, int *y2, int *xy, PlanarFrame &msk, 
							  PlanarFrame &dst, const int pitch)
{
	unsigned char *mskp = msk.GetPtr(plane);
	const int msk_pitch = msk.GetPitch(plane);
	const int height = msk.GetHeight(plane);
	const int width = msk.GetWidth(plane);
	unsigned char *dstp = dst.GetPtr(plane);
	const int dst_pitch = dst.GetPitch(plane);
	mskp += (8-field)*msk_pitch;
	dstp += (8-field)*dst_pitch;
	unsigned char *dstpp = dstp-dst_pitch;
	unsigned char *dstpn = dstp+dst_pitch;
	x2 += pitch*3;
	y2 += pitch*3;
	xy += pitch*3;
	int *x2n = x2+pitch;
	int *y2n = y2+pitch;
	int *xyn = xy+pitch;
	for (int y=8-field; y<height-7; y+=2)
	{
		for (int x=4; x<width-4; ++x)
		{
			if (mskp[x] == 255 || mskp[x] == 128) continue;
			const int c1 = int(x2[x]*y2[x]-xy[x]*xy[x]-0.09*(x2[x]+y2[x])*(x2[x]+y2[x]));
			const int c2 = int(x2n[x]*y2n[x]-xyn[x]*xyn[x]-0.09*(x2n[x]+y2n[x])*(x2n[x]+y2n[x]));
			if (c1 > 775 || c2 > 775)
				dstp[x] = (dstpp[x]+dstpn[x]+1)>>1;
		}
		mskp += msk_pitch*2;
		dstpp += dst_pitch*2;
		dstp += dst_pitch*2;
		dstpn += dst_pitch*2;
		x2 += pitch;
		x2n += pitch;
		y2 += pitch;
		y2n += pitch;
		xy += pitch;
		xyn += pitch;
	}
}

void EEDI2::gaussianBlur1(const int plane, PlanarFrame &src, PlanarFrame &tmp, PlanarFrame &dst)
{
	unsigned char *srcp = src.GetPtr(plane);
	const int src_pitcht = src.GetPitch(plane);
	const int width = src.GetWidth(plane);
	const int height = src.GetHeight(plane);
	unsigned char *dstp = tmp.GetPtr(plane);
	const int tmp_pitch = tmp.GetPitch(plane);
	for (int y=0; y<height; ++y)
	{
		dstp[0] = (srcp[3]*582 + srcp[2]*7078 + srcp[1]*31724 + 
				srcp[0]*26152 + 32768)>>16;
		dstp[1] = (srcp[4]*582 + srcp[3]*7078 +
				(srcp[0]+srcp[2])*15862 + srcp[1]*26152 + 32768)>>16;
		dstp[2] = (srcp[5]*582 + (srcp[0]+srcp[4])*3539 +
				(srcp[1]+srcp[3])*15862 + srcp[2]*26152 + 32768)>>16;
		for (int x=3; x<width-3; ++x)
		{
			dstp[x] = ((srcp[x-3]+srcp[x+3])*291 + (srcp[x-2]+srcp[x+2])*3539 +
					(srcp[x-1]+srcp[x+1])*15862 + srcp[x]*26152 + 32768)>>16;
		}
		dstp[x] = (srcp[x-3]*582 + (srcp[x-2]+srcp[x+2])*3539 +
				(srcp[x-1]+srcp[x+1])*15862 + srcp[x]*26152 + 32768)>>16; ++x;
		dstp[x] = (srcp[x-3]*582 + srcp[x-2]*7078 +
				(srcp[x-1]+srcp[x+1])*15862 + srcp[x]*26152 + 32768)>>16; ++x;
		dstp[x] = (srcp[x-3]*582 + srcp[x-2]*7078 + srcp[x-1]*31724 + 
				srcp[x]*26152 + 32768)>>16;
		srcp += src_pitcht;
		dstp += tmp_pitch;
	}
	srcp = tmp.GetPtr(plane);
	dstp = dst.GetPtr(plane);
	const int src_pitch = tmp.GetPitch(plane);
	const int dst_pitch = dst.GetPitch(plane);
	unsigned char *src3p = srcp-src_pitch*3;
	unsigned char *src2p = srcp-src_pitch*2;
	unsigned char *srcpp = srcp-src_pitch;
	unsigned char *srcpn = srcp+src_pitch;
	unsigned char *src2n = srcp+src_pitch*2;
	unsigned char *src3n = srcp+src_pitch*3;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3n[x]*582 + src2n[x]*7078 + srcpn[x]*31724 + 
			srcp[x]*26152 + 32768)>>16;
	}
	src3p += src_pitch;
	src2p += src_pitch;
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	src2n += src_pitch;
	src3n += src_pitch;
	dstp += dst_pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3n[x]*582 + src2n[x]*7078 + (srcpp[x]+srcpn[x])*15862 +
				srcp[x]*26152 + 32768)>>16;
	}
	src3p += src_pitch;
	src2p += src_pitch;
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	src2n += src_pitch;
	src3n += src_pitch;
	dstp += dst_pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3n[x]*582 + (src2p[x]+src2n[x])*3539 + 
				(srcpp[x]+srcpn[x])*15862 + srcp[x]*26152 + 32768)>>16;
	}
	src3p += src_pitch;
	src2p += src_pitch;
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	src2n += src_pitch;
	src3n += src_pitch;
	dstp += dst_pitch;
	for (int y=3; y<height-3; ++y)
	{
		for (int x=0; x<width; ++x)
		{
			dstp[x] = ((src3p[x]+src3n[x])*291 +
					(src2p[x]+src2n[x])*3539 + (srcpp[x]+srcpn[x])*15862 +
					srcp[x]*26152 + 32768)>>16;
		}
		src3p += src_pitch;
		src2p += src_pitch;
		srcpp += src_pitch;
		srcp += src_pitch;
		srcpn += src_pitch;
		src2n += src_pitch;
		src3n += src_pitch;
		dstp += dst_pitch;
	}
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3p[x]*582 + (src2p[x]+src2n[x])*3539 + (srcpp[x]+srcpn[x])*15862 +
				srcp[x]*26152 + 32768)>>16;
	}
	src3p += src_pitch;
	src2p += src_pitch;
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	src2n += src_pitch;
	src3n += src_pitch;
	dstp += dst_pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3p[x]*582 + src2p[x]*7078 + (srcpp[x]+srcpn[x])*15862 +
				srcp[x]*26152 + 32768)>>16;
	}
	src3p += src_pitch;
	src2p += src_pitch;
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	src2n += src_pitch;
	src3n += src_pitch;
	dstp += dst_pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src3p[x]*582 + src2p[x]*7078 + srcpp[x]*31724 + 
			srcp[x]*26152 + 32768)>>16;
	}
}

void EEDI2::gaussianBlurSqrt2(int *src, int *tmp, int *dst, const int width, const int height, 
		const int pitch)
{
	int *srcp = src;
	int *dstp = tmp;
	for (int y=0; y<height; ++y)
	{
		int x=0;
		dstp[x] = (srcp[x+4]*678 + srcp[x+3]*3902 + srcp[x+2]*13618 + srcp[x+1]*28830 + 
				srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x+4]*678 + srcp[x+3]*3902 + srcp[x+2]*13618 + 
			(srcp[x-1]+srcp[x+1])*14415 + srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x+4]*678 + srcp[x+3]*3902 + 
			(srcp[x-2]+srcp[x+2])*6809 + (srcp[x-1]+srcp[x+1])*14415 + 
			srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x+4]*678 + (srcp[x-3]+srcp[x+3])*1951 + 
			(srcp[x-2]+srcp[x+2])*6809 + (srcp[x-1]+srcp[x+1])*14415 + 
			srcp[x]*18508 + 32768)>>16;
		for (x=4; x<width-4; ++x)
		{
			dstp[x] = ((srcp[x-4]+srcp[x+4])*339 + (srcp[x-3]+srcp[x+3])*1951 + 
					(srcp[x-2]+srcp[x+2])*6809 + (srcp[x-1]+srcp[x+1])*14415 + 
					srcp[x]*18508 + 32768)>>16;
		}
		dstp[x] = (srcp[x-4]*678 + (srcp[x-3]+srcp[x+3])*1951 + 
				(srcp[x-2]+srcp[x+2])*6809 + (srcp[x-1]+srcp[x+1])*14415 + 
				srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x-4]*678 + srcp[x-3]*3902 + 
			(srcp[x-2]+srcp[x+2])*6809 + (srcp[x-1]+srcp[x+1])*14415 + 
			srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x-4]*678 + srcp[x+3]*3902 + srcp[x-2]*13618 + 
			(srcp[x-1]+srcp[x+1])*14415 + srcp[x]*18508 + 32768)>>16; ++x;
		dstp[x] = (srcp[x-4]*678 + srcp[x-3]*3902 + srcp[x-2]*13618 + 
			srcp[x-1]*28830 + srcp[x]*18508 + 32768)>>16;
		srcp += pitch;
		dstp += pitch;
	}
	dstp = dst;
	srcp = tmp;
	int *src4p = srcp-pitch*4;
	int *src3p = srcp-pitch*3;
	int *src2p = srcp-pitch*2;
	int *srcpp = srcp-pitch;
	int *srcpn = srcp+pitch;
	int *src2n = srcp+pitch*2;
	int *src3n = srcp+pitch*3;
	int *src4n = srcp+pitch*4;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4n[x]*678 + src3n[x]*3902 + src2n[x]*13618 + srcpn[x]*28830 +
				srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4n[x]*678 + src3n[x]*3902 + src2n[x]*13618 + 
			(srcpp[x]+srcpn[x])*14415 + srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4n[x]*678 + src3n[x]*3902 + (src2p[x]+src2n[x])*6809 + 
			(srcpp[x]+srcpn[x])*14415 + srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4n[x]*678 + (src3p[x]+src3n[x])*1951 +
				(src2p[x]+src2n[x])*6809 + (srcpp[x]+srcpn[x])*14415 +
				srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int y=4; y<height-4; ++y)
	{
		for (int x=0; x<width; ++x)
		{
			dstp[x] = ((src4p[x]+src4n[x])*339 + (src3p[x]+src3n[x])*1951 +
					(src2p[x]+src2n[x])*6809 + (srcpp[x]+srcpn[x])*14415 +
					srcp[x]*18508 + 32768)>>18;
		}
		src4p += pitch;
		src3p += pitch;
		src2p += pitch;
		srcpp += pitch;
		srcp += pitch;
		srcpn += pitch;
		src2n += pitch;
		src3n += pitch;
		src4n += pitch;
		dstp += pitch;
	}
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4p[x]*678 + (src3p[x]+src3n[x])*1951 +
				(src2p[x]+src2n[x])*6809 + (srcpp[x]+srcpn[x])*14415 +
				srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4p[x]*678 + src3p[x]*3902 +
				(src2p[x]+src2n[x])*6809 + (srcpp[x]+srcpn[x])*14415 +
				srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4p[x]*678 + src3p[x]*3902 +
				src2p[x]*13618 + (srcpp[x]+srcpn[x])*14415 +
				srcp[x]*18508 + 32768)>>18;
	}
	src4p += pitch;
	src3p += pitch;
	src2p += pitch;
	srcpp += pitch;
	srcp += pitch;
	srcpn += pitch;
	src2n += pitch;
	src3n += pitch;
	src4n += pitch;
	dstp += pitch;
	for (int x=0; x<width; ++x)
	{
		dstp[x] = (src4p[x]*678 + src3p[x]*3902 +
				src2p[x]*13618 + srcpp[x]*28830 +
				srcp[x]*18508 + 32768)>>18;
	}
}

void EEDI2::calcDerivatives(const int plane, PlanarFrame &src, int *x2, int *y2, int *xy)
{
	unsigned char *srcp = src.GetPtr(plane);
	const int src_pitch = src.GetPitch(plane);
	const int width = src.GetWidth(plane);
	const int height = src.GetHeight(plane);
	unsigned char *srcpp = srcp-src_pitch;
	unsigned char *srcpn = srcp+src_pitch;
	{
		const int Ix = srcp[1]-srcp[0];
		const int Iy = srcp[0]-srcpn[0];
		x2[0] = (Ix*Ix)>>1;
		y2[0] = (Iy*Iy)>>1;
		xy[0] = (Ix*Iy)>>1;
	}
	for (int x=1; x<width-1; ++x)
	{
		const int Ix = srcp[x+1]-srcp[x-1];
		const int Iy = srcp[x]-srcpn[x];
		x2[x] = (Ix*Ix)>>1;
		y2[x] = (Iy*Iy)>>1;
		xy[x] = (Ix*Iy)>>1;
	}
	{
		const int Ix = srcp[x]-srcp[x-1];
		const int Iy = srcp[x]-srcpn[x];
		x2[x] = (Ix*Ix)>>1;
		y2[x] = (Iy*Iy)>>1;
		xy[x] = (Ix*Iy)>>1;
	}
	srcpp += src_pitch;
	srcp += src_pitch;
	srcpn += src_pitch;
	x2 += src_pitch;
	y2 += src_pitch;
	xy += src_pitch;
	for (int y=1; y<height-1; ++y)
	{
		{
			const int Ix = srcp[1]-srcp[0];
			const int Iy = srcpp[0]-srcpn[0];
			x2[0] = (Ix*Ix)>>1;
			y2[0] = (Iy*Iy)>>1;
			xy[0] = (Ix*Iy)>>1;
		}
		for (int x=1; x<width-1; ++x)
		{
			const int Ix = srcp[x+1]-srcp[x-1];
			const int Iy = srcpp[x]-srcpn[x];
			x2[x] = (Ix*Ix)>>1;
			y2[x] = (Iy*Iy)>>1;
			xy[x] = (Ix*Iy)>>1;
		}
		{
			const int Ix = srcp[x]-srcp[x-1];
			const int Iy = srcpp[x]-srcpn[x];
			x2[x] = (Ix*Ix)>>1;
			y2[x] = (Iy*Iy)>>1;
			xy[x] = (Ix*Iy)>>1;
		}
		srcpp += src_pitch;
		srcp += src_pitch;
		srcpn += src_pitch;
		x2 += src_pitch;
		y2 += src_pitch;
		xy += src_pitch;
	}
	{
		const int Ix = srcp[1]-srcp[0];
		const int Iy = srcpp[0]-srcp[0];
		x2[0] = (Ix*Ix)>>1;
		y2[0] = (Iy*Iy)>>1;
		xy[0] = (Ix*Iy)>>1;
	}
	for (int x=1; x<width-1; ++x)
	{
		const int Ix = srcp[x+1]-srcp[x-1];
		const int Iy = srcpp[x]-srcp[x];
		x2[x] = (Ix*Ix)>>1;
		y2[x] = (Iy*Iy)>>1;
		xy[x] = (Ix*Iy)>>1;
	}
	{
		const int Ix = srcp[x]-srcp[x-1];
		const int Iy = srcpp[x]-srcp[x];
		x2[x] = (Ix*Ix)>>1;
		y2[x] = (Iy*Iy)>>1;
		xy[x] = (Ix*Iy)>>1;
	}
}

AVSValue __cdecl Create_EEDI2(AVSValue args, void* user_data, IScriptEnvironment* env) 
{
	return new EEDI2(args[0].AsClip(),args[1].AsInt(10),args[2].AsInt(20),
		args[3].AsInt(20),args[4].AsInt(2),args[5].AsInt(4),args[6].AsInt(24),
		args[7].AsInt(-1),args[8].AsInt(0),args[9].AsInt(50),args[10].AsInt(1),
		env);
}

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit2(IScriptEnvironment* env) 
{
    env->AddFunction("EEDI2", "c[mthresh]i[lthresh]i[vthresh]i[estr]i[dstr]i[maxd]i[field]i" \
		"[map]i[nt]i[pp]i", Create_EEDI2, 0);
    return 0;
}