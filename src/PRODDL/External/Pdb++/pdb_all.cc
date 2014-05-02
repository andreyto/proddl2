//
//	Copyright (c) 1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb++.cc,v 1.6 1994/12/13 22:41:17 gregc Exp $
//

#include "PRODDL/External/Pdb++/pdb++.hpp"

#include <cstring>

namespace PDBPP {

using namespace std;

void
PDB::type(RecordType t)
{
	if (t == UNKNOWN) {
		// optimize default case (skip memset())
		rType = t;
		unknown.junk[0] = '\0';
		return;
	}
	memset(this, 0, sizeof *this);
	rType = t;
	switch (t) {
	default:
		break;
	case ATOM:
		atom.occupancy = 1.0;
		break;
	}
}

int
PDB::byteCmp(const PDB &l, const PDB &r)
{
	return memcmp(&l, &r, sizeof (PDB));
}

} // namespace PDBPP
//
//	Copyright (c) 1989,1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_chars.cc,v 1.8 95/02/28 14:09:52 gregc Exp $
//
//	subroutine for writing PDB format files
//


# include	<ctype.h>

extern "C" int	sprintf(char *, const char *, ...);

namespace PDBPP {

using namespace std;

static char const * const pdbRecordFormatW[PDB::NUM_TYPES] = {
#include "PRODDL/External/Pdb++/write_format.i"
};

static char const * const pdbrun5W[PDB::NUM_USER] = {
#include "PRODDL/External/Pdb++/pdbrun5_write.i"
};

static char const * const pdbrun6W[PDB::NUM_USER] = {
#include "PRODDL/External/Pdb++/pdbrun6_write.i"
};

const char *
PDB::chars(void) const
{
	static char	buf[BufLen];
	const char	*fmt;
	const Sheet	*sh;
	const Residue	*shr0, *shr1, *sha0, *sha1;
	int		count;

	// convert C structure to pdb record

	if (rType < USER_PDBRUN)
		fmt = pdbRecordFormatW[rType];
	else if (pdbrunOutputVersion < 6)
		fmt = pdbrun5W[rType - USER_PDBRUN];
	else
		fmt = pdbrun6W[rType - USER_PDBRUN];
	switch (rType) {

	case UNKNOWN:
		count = sprintf(buf, fmt, unknown.junk);
		break;

	case AGGRGT:
		count = sprintf(buf, fmt, aggrgt.serialNum,
			aggrgt.numComponents, aggrgt.cmpontSerialNums[0],
			aggrgt.cmpontSerialNums[1],
			aggrgt.cmpontSerialNums[2],
			aggrgt.cmpontSerialNums[3],
			aggrgt.cmpontSerialNums[4],
			aggrgt.cmpontSerialNums[5],
			aggrgt.cmpontSerialNums[6],
			aggrgt.cmpontSerialNums[7],
			aggrgt.cmpontSerialNums[8],
			aggrgt.cmpontSerialNums[9],
			aggrgt.cmpontSerialNums[10],
			aggrgt.cmpontSerialNums[11],
			aggrgt.cmpontSerialNums[12],
			aggrgt.cmpontSerialNums[13]);
		break;

	case AGRDES:
	case CMPDES:
	case FTNOTE:
	case MTXDES:
	case REMARK:
	case SYMDES:
		count = sprintf(buf, fmt, agrdes.num, agrdes.text);
		break;

	case ANISOU:
	case SIGUIJ:
		count = sprintf(buf, fmt, anisou.serialNum, anisou.name,
			anisou.altLoc, anisou.residue.name,
			anisou.residue.chainId, anisou.residue.seqNum,
			anisou.residue.insertCode, anisou.u[0], anisou.u[1],
			anisou.u[2], anisou.u[3], anisou.u[4], anisou.u[5]);
		break;

	case ATOM:
	case HETATM:
	case SIGATM:
		count = sprintf(buf, fmt, atom.serialNum, atom.name,
			atom.altLoc, atom.residue.name, atom.residue.chainId,
			atom.residue.seqNum, atom.residue.insertCode,
			atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.occupancy,
			atom.tempFactor, atom.ftnoteNum);
		break;

	case AUTHOR:
	case COMPND:
	case JRNL:
	case SOURCE:
	case EXPDTA:
		count = sprintf(buf, fmt, author.continuation, author.data);
		break;

	case CONECT:
		count = sprintf(buf, fmt, conect.serialNum,
			conect.covalent[0], conect.covalent[1],
			conect.covalent[2], conect.covalent[3],
			conect.bonds[0].hydrogen[0],
			conect.bonds[0].hydrogen[1], conect.bonds[0].salt,
			conect.bonds[1].hydrogen[0],
			conect.bonds[1].hydrogen[1], conect.bonds[1].salt);
		break;

	case CMPONT:
		count = sprintf(buf, fmt, cmpont.seqNum,
			cmpont.residues[0].name, cmpont.residues[0].chainId,
			cmpont.residues[0].seqNum,
			cmpont.residues[0].insertCode,
			cmpont.residues[1].name, cmpont.residues[1].chainId,
			cmpont.residues[1].seqNum,
			cmpont.residues[1].insertCode);
		break;

	case CRYST1:
		count = sprintf(buf, fmt, cryst1.a, cryst1.b, cryst1.c,
			cryst1.alpha, cryst1.beta, cryst1.gamma,
			cryst1.spaceGrp, cryst1.z);
		break;

	case END:
	case ENDMDL:
		count = sprintf(buf, fmt);
		break;

	case FORMUL:
		count = sprintf(buf, fmt, formul.component, formul.hetId,
			formul.continuation, formul.exclude, formul.formula);
		break;

	case HEADER:
		count = sprintf(buf, fmt, header.classification,
			header.timestamp, header.type, header.id);
		break;

	case HELIX:
		count = sprintf(buf, fmt, helix.serialNum, helix.id,
			helix.residues[0].name, helix.residues[0].chainId,
			helix.residues[0].seqNum,
			helix.residues[0].insertCode, helix.residues[1].name,
			helix.residues[1].chainId, helix.residues[1].seqNum,
			helix.residues[1].insertCode, helix.type,
			helix.comment);
		break;

	case HET:
		count = sprintf(buf, fmt, het.hetGrp.name,
			het.hetGrp.chainId, het.hetGrp.seqNum,
			het.hetGrp.insertCode, het.numAtoms, het.text);
		break;

	case MASTER:
		count = sprintf(buf, fmt, master.numRemark, master.numFtnote,
			master.numHet, master.numHelix, master.numSheet,
			master.numTurn, master.numSite, master.numTransform,
			master.numCoordinate, master.numTer,
			master.numConect, master.numSeqres);
		break;

	case MODEL:
		count = sprintf(buf, fmt, model.num);
		break;

	case MTRIX:
		count = sprintf(buf, fmt, mtrix.rowNum, mtrix.serialNum,
			mtrix.m1, mtrix.m2, mtrix.m3, mtrix.v, mtrix.given);
		break;

	case OBSLTE:
		count = sprintf(buf, fmt, obslte.continuation, obslte.timestamp,
			obslte.oldId, obslte.idMap[0], obslte.idMap[1],
			obslte.idMap[2], obslte.idMap[3], obslte.idMap[4],
			obslte.idMap[2], obslte.idMap[6], obslte.idMap[7]);
		break;

	case ORIGX:
		count = sprintf(buf, fmt, origx.rowNum, origx.o1, origx.o2,
			origx.o3, origx.t);
		break;

	case REVDAT:
		count = sprintf(buf, fmt, revdat.modification,
			revdat.continuation, revdat.timestamp, revdat.id,
			revdat.modType, revdat.corrections);
		break;

	case SCALE:
		count = sprintf(buf, fmt, scale.rowNum, scale.s1, scale.s2,
			scale.s3, scale.u);
		break;

	case SEQRES:
		count = sprintf(buf, fmt, seqres.serialNum, seqres.chainId,
			seqres.count, seqres.names[0], seqres.names[1],
			seqres.names[2], seqres.names[3], seqres.names[4],
			seqres.names[5], seqres.names[6], seqres.names[7],
			seqres.names[8], seqres.names[9], seqres.names[10],
			seqres.names[11], seqres.names[12]);
		break;

	case SHEET:
		sh = &sheet;
		shr0 = &sh->residues[0];
		shr1 = &sh->residues[1];
		sha0 = &sh->atoms[0].residue;
		sha1 = &sh->atoms[1].residue;
		count = sprintf(buf, fmt, sh->strandNum, sh->id, sh->count,
			shr0->name, shr0->chainId, shr0->seqNum,
			shr0->insertCode, shr1->name, shr1->chainId,
			shr1->seqNum, shr1->insertCode, sh->sense,
			sh->atoms[0].name, sha0->name, sha0->chainId,
			sha0->seqNum, sha0->insertCode, sh->atoms[1].name,
			sha1->name, sha1->chainId, sha1->seqNum,
			sha1->insertCode);
		break;

	case SITE:
		shr0 = &site.residues[0];
		shr1 = &site.residues[1];
		sha0 = &site.residues[2];
		sha1 = &site.residues[3];
		count = sprintf(buf, fmt, site.seqNum, site.id, site.count,
			shr0->name, shr0->chainId, shr0->seqNum,
			shr0->insertCode,
			shr1->name, shr1->chainId, shr1->seqNum,
			shr1->insertCode,
			sha0->name, sha0->chainId, sha0->seqNum,
			sha0->insertCode,
			sha1->name, sha1->chainId, sha1->seqNum,
			sha1->insertCode);
		break;

	case SPRSDE:
		count = sprintf(buf, fmt, sprsde.continuation, sprsde.timestamp,
			sprsde.id, sprsde.supersede[0], sprsde.supersede[1],
			sprsde.supersede[2], sprsde.supersede[3],
			sprsde.supersede[4], sprsde.supersede[5],
			sprsde.supersede[6], sprsde.supersede[7]);
		break;

	case SSBOND:
		count = sprintf(buf, fmt, ssbond.seqNum,
			ssbond.residues[0].name, ssbond.residues[0].chainId,
			ssbond.residues[0].seqNum,
			ssbond.residues[0].insertCode,
			ssbond.residues[1].name, ssbond.residues[1].chainId,
			ssbond.residues[1].seqNum,
			ssbond.residues[1].insertCode, ssbond.comment);
		break;

	case SYMOP:
		count = sprintf(buf, fmt, symop.rowNum, symop.serialNum,
			symop.s1, symop.s2, symop.s3, symop.t);
		break;

	case TER:
		count = sprintf(buf, fmt, ter.serialNum, ter.residue.name,
			ter.residue.chainId, ter.residue.seqNum,
			ter.residue.insertCode);
		break;

	case TRNSFM:
		count = sprintf(buf, fmt, trnsfm.resultSerialNum,
			trnsfm.applySerialNum, trnsfm.sourceSerialNum);
		break;

	case TURN:
		count = sprintf(buf, fmt, turn.seqNum, turn.id,
			turn.residues[0].name, turn.residues[0].chainId,
			turn.residues[0].seqNum, turn.residues[0].insertCode,
			turn.residues[1].name, turn.residues[1].chainId,
			turn.residues[1].seqNum, turn.residues[1].insertCode,
			turn.comment);
		break;

	case TVECT:
		count = sprintf(buf, fmt, tvect.serialNum, tvect.t1, tvect.t2,
			tvect.t3, tvect.comment);
		break;

	case USER:
		count = sprintf(buf, fmt, user.subtype, user.text);
		break;

	case USER_PDBRUN:
		count = sprintf(buf, fmt, userPdbrun.version);
		pdbrunOutputVersion = userPdbrun.version;
		break;

	case USER_EYEPOS:
		count = sprintf(buf, fmt, userEyePos.xyz[0], userEyePos.xyz[1],
			userEyePos.xyz[2]);
		break;

	case USER_ATPOS:
		count = sprintf(buf, fmt, userAtPos.xyz[0], userAtPos.xyz[1],
			userAtPos.xyz[2]);
		break;

	case USER_WINDOW:
		count = sprintf(buf, fmt, userWindow.left, userWindow.right,
			userWindow.bottom, userWindow.top, userWindow.hither,
			userWindow.yon);
		break;

	case USER_FOCUS:
		count = sprintf(buf, fmt, userFocus.focus);
		break;

	case USER_VIEWPORT:
		count = sprintf(buf, fmt, userViewport.xmin, userViewport.xmax,
			userViewport.ymin, userViewport.ymax);
		break;

	case USER_BGCOLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userBgColor.rgb[0],
				userBgColor.rgb[1], userBgColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userBgColor.rgb[0],
				userBgColor.rgb[1], userBgColor.rgb[2]);
		break;

	case USER_ANGLE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userAngle.which,
				userAngle.atom0, userAngle.atom1,
				userAngle.atom2, userAngle.atom3,
				userAngle.angle);
		else
			count = sprintf(buf, fmt, userAngle.atom0,
				userAngle.atom1, userAngle.atom2,
				userAngle.atom3, userAngle.angle);
		break;

	case USER_DISTANCE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userDistance.which,
				userDistance.atom0, userDistance.atom1,
				userDistance.distance);
		else
			count = sprintf(buf, fmt, userDistance.atom0,
				userDistance.atom1, userDistance.distance);
		break;

	case USER_FILE:
		if (pdbrunOutputVersion < 6)
			count = sprintf(buf, fmt, userFile.filename);
		else
			count = sprintf(buf, fmt, userFile.model,
							userFile.filename);
		break;

	case USER_MARKNAME:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userMarkname.markname);
		break;

	case USER_MARK:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userMark.markname);
		break;

	case USER_CNAME:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userCName.name,
				userCName.rgb[0], userCName.rgb[1],
				userCName.rgb[2]);
		else
			count = sprintf(buf, fmt, userCName.rgb[0],
				userCName.rgb[1], userCName.rgb[2],
				userCName.name);
		break;

	case USER_COLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userColor.spec,
				userColor.rgb[0], userColor.rgb[1],
				userColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userColor.rgb[0],
				userColor.rgb[1], userColor.rgb[2],
				userColor.spec);
		break;

	case USER_RADIUS:
		count = sprintf(buf, fmt, userRadius.radius);
		break;

	case USER_OBJECT:
		count = sprintf(buf, fmt, userObject.model);
		break;

	case USER_ENDOBJ:
		count = sprintf(buf, fmt, userEndObj.model);
		break;

	case USER_CHAIN:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userChain.atom0,
				userChain.atom1);
		else
			count = sprintf(buf, fmt, userChain.atom0,
				userChain.atom1);
		break;

	case USER_GFX_BEGIN:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else if (userGfxBegin.primitive == GFX_UNKNOWN)
			count = sprintf(buf, fmt, userGfxBegin.unknown);
		else
			count = sprintf(buf, fmt,
					gfxChars(userGfxBegin.primitive));
		break;

	case USER_GFX_END:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt);
		break;

	case USER_GFX_COLOR:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxColor.spec,
				userGfxColor.rgb[0], userGfxColor.rgb[1],
				userGfxColor.rgb[2]);
		else
			count = sprintf(buf, fmt, userGfxColor.rgb[0],
				userGfxColor.rgb[1], userGfxColor.rgb[2],
				userGfxColor.spec);
		break;

	case USER_GFX_NORMAL:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxNormal.xyz[0],
				userGfxNormal.xyz[1], userGfxNormal.xyz[2]);
		break;

	case USER_GFX_VERTEX:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxVertex.xyz[0],
				userGfxVertex.xyz[1], userGfxVertex.xyz[2]);
		break;

	case USER_GFX_FONT:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxFont.name,
				userGfxFont.size);
		else
			count = sprintf(buf, fmt, userGfxFont.size,
				userGfxFont.name);
		break;

	case USER_GFX_TEXTPOS:
		if (pdbrunOutputVersion < 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxTextPos.xyz[0],
				userGfxTextPos.xyz[1], userGfxTextPos.xyz[2]);
		break;

	case USER_GFX_LABEL:
		if (pdbrunOutputVersion < 6)
			count = ::sprintf(buf, fmt, userGfxLabel.xyz[0],
				userGfxLabel.xyz[1], userGfxLabel.xyz[2],
				userGfxLabel.text);
		else
			count = sprintf(buf, fmt, userGfxLabel.text);
		break;

	case USER_GFX_MOVE:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxMove.xyz[0],
				userGfxMove.xyz[1], userGfxMove.xyz[2]);
		break;

	case USER_GFX_DRAW:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxDraw.xyz[0],
				userGfxDraw.xyz[1], userGfxDraw.xyz[2]);
		break;

	case USER_GFX_MARKER:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxMarker.xyz[0],
				userGfxMarker.xyz[1], userGfxMarker.xyz[2]);
		break;

	case USER_GFX_POINT:
		if (pdbrunOutputVersion >= 6)
			count = 0;
		else
			count = sprintf(buf, fmt, userGfxPoint.xyz[0],
				userGfxPoint.xyz[1], userGfxPoint.xyz[2]);
		break;

	default:
		count = sprintf(buf, "unknown pdb record #%d", rType);
		break;
	}

	// find last non-blank in buf, and shorten it
	while (count > 1 && isspace(buf[count - 1]))
		count -= 1;
	buf[count] = '\0';
	return buf;
}

} // namespace PDBPP
//
//	Copyright (c) 1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdbinput.cc,v 1.1 94/02/22 17:09:27 gregc Exp $
//


namespace std {

istream &
operator>>(istream &s, PDBPP::PDB &p)
{
    using namespace PDBPP;
	char	buf[4 * PDB::BufLen];

	s.getline(buf, 4 * PDB::BufLen);
	if (s)
		p = PDB(buf);
	return s;
}

} // namespace std
//
//	Copyright (c) 1989, 1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_read.cc,v 1.8 95/02/28 14:09:28 gregc Exp $
//
//	subroutine for reading PDB format files
//


# include	<string.h>

extern "C" int	sscanf(const char *, const char *, ...);

//
//	for each pdb record type there is a format reading in the
//	record values and for printing them out.
//
//	The actual format of a line written, is the print format
//	followed by blank padding to 72 characters, followed by
//	8 characters of file and line information.
//

namespace PDBPP {

using namespace std;

static char const * const pdbRecordFormatR[PDB::NUM_TYPES] = {
#include "PRODDL/External/Pdb++/read_format.i"
};

static char const * const pdbrun5R[PDB::NUM_USER] = {
#include "PRODDL/External/Pdb++/pdbrun5_read.i"
};

static char const * const pdbrun6R[PDB::NUM_USER] = {
#include "PRODDL/External/Pdb++/pdbrun6_read.i"
};

PDB::PDB(const char *buf)
{
	const char	*fmt;
	Sheet		*sh;
	Residue		*sha0, *sha1;

	// convert pdb record to C structure

	memset(this, 0, sizeof *this);
	rType = getType(buf);
	if (rType < USER_PDBRUN)
		fmt = pdbRecordFormatR[rType];
	else if (pdbrunInputVersion < 6)
		fmt = pdbrun5R[rType - USER_PDBRUN];
	else
		fmt = pdbrun6R[rType - USER_PDBRUN];
	switch (rType) {

	default:
	case UNKNOWN:
unknown:
		rType = UNKNOWN;		// in case of goto
		(void) sprintf(unknown.junk, "%72s", buf);
		break;

	case AGGRGT:
		if (0 > sscanf(buf, fmt, &aggrgt.serialNum,
				&aggrgt.numComponents,
				&aggrgt.cmpontSerialNums[0],
				&aggrgt.cmpontSerialNums[1],
				&aggrgt.cmpontSerialNums[2],
				&aggrgt.cmpontSerialNums[3],
				&aggrgt.cmpontSerialNums[4],
				&aggrgt.cmpontSerialNums[5],
				&aggrgt.cmpontSerialNums[6],
				&aggrgt.cmpontSerialNums[7],
				&aggrgt.cmpontSerialNums[8],
				&aggrgt.cmpontSerialNums[9],
				&aggrgt.cmpontSerialNums[10],
				&aggrgt.cmpontSerialNums[11],
				&aggrgt.cmpontSerialNums[12],
				&aggrgt.cmpontSerialNums[13]))
			goto unknown;
		break;

	case AGRDES:
	case CMPDES:
	case FTNOTE:
	case MTXDES:
	case REMARK:
	case SYMDES:
		if (0 > sscanf(buf, fmt, &agrdes.num, agrdes.text))
			goto unknown;
		break;

	case ANISOU:
	case SIGUIJ:
		if (0 > sscanf(buf, fmt, &anisou.serialNum, anisou.name,
				&anisou.altLoc, anisou.residue.name,
				&anisou.residue.chainId,
				&anisou.residue.seqNum,
				&anisou.residue.insertCode,
				&anisou.u[0], &anisou.u[1], &anisou.u[2],
				&anisou.u[3], &anisou.u[4], &anisou.u[5]))
			goto unknown;
		break;

	case ATOM:
	case HETATM:
	case SIGATM:
	    //DEBUG:
	    //std::cout << "sscanf ATOM\n";
		if (0 > sscanf(buf, fmt, &atom.serialNum, atom.name,
				&atom.altLoc, atom.residue.name,
				&atom.residue.chainId, &atom.residue.seqNum,
				&atom.residue.insertCode, &atom.xyz[0],
				&atom.xyz[1], &atom.xyz[2], &atom.occupancy)) {
				//,
				//&atom.tempFactor, &atom.ftnoteNum)) {
			//std::cout << "goto unknown\n";
			goto unknown;
			}
		break;

	case AUTHOR:
	case COMPND:
	case EXPDTA:
	case JRNL:
	case SOURCE:
		if (0 > sscanf(buf, fmt, &author.continuation, author.data))
			goto unknown;
		break;

	case CONECT:
		if (0 > sscanf(buf, fmt, &conect.serialNum,
				&conect.covalent[0], &conect.covalent[1],
				&conect.covalent[2], &conect.covalent[3],
				&conect.bonds[0].hydrogen[0],
				&conect.bonds[0].hydrogen[1],
				&conect.bonds[0].salt,
				&conect.bonds[1].hydrogen[0],
				&conect.bonds[1].hydrogen[1],
				&conect.bonds[1].salt))
			goto unknown;
		break;

	case CMPONT:
		if (0 > sscanf(buf, fmt, &cmpont.seqNum,
				cmpont.residues[0].name,
				&cmpont.residues[0].chainId,
				&cmpont.residues[0].seqNum,
				&cmpont.residues[0].insertCode,
				cmpont.residues[1].name,
				&cmpont.residues[1].chainId,
				&cmpont.residues[1].seqNum,
				&cmpont.residues[1].insertCode))
			goto unknown;
		break;

	case CRYST1:
		if (0 > sscanf(buf, fmt, &cryst1.a, &cryst1.b, &cryst1.c,
				&cryst1.alpha, &cryst1.beta, &cryst1.gamma,
				cryst1.spaceGrp, &cryst1.z))
			goto unknown;
		break;

	case END:
	case ENDMDL:
		break;

	case FORMUL:
		if (0 > sscanf(buf, fmt, &formul.component, formul.hetId,
				&formul.continuation, &formul.exclude,
				formul.formula))
			goto unknown;
		break;

	case HEADER:
		if (0 > sscanf(buf, fmt, header.classification,
				header.timestamp, &header.type, header.id))
			goto unknown;
		break;

	case HELIX:
		if (0 > sscanf(buf, fmt, &helix.serialNum, helix.id,
				helix.residues[0].name,
				&helix.residues[0].chainId,
				&helix.residues[0].seqNum,
				&helix.residues[0].insertCode,
				helix.residues[1].name,
				&helix.residues[1].chainId,
				&helix.residues[1].seqNum,
				&helix.residues[1].insertCode,
				&helix.type, helix.comment))
			goto unknown;
		break;

	case HET:
		if (0 > sscanf(buf, fmt, het.hetGrp.name,
				&het.hetGrp.chainId, &het.hetGrp.seqNum,
				&het.hetGrp.insertCode, &het.numAtoms,
				het.text))
			goto unknown;
		break;

	case MASTER:
		if (0 > sscanf(buf, fmt, &master.numRemark, &master.numFtnote,
				&master.numHet, &master.numHelix,
				&master.numSheet, &master.numTurn,
				&master.numSite, &master.numTransform,
				&master.numCoordinate, &master.numTer,
				&master.numConect, &master.numSeqres))
			goto unknown;
		break;

	case MODEL:
		if (0 > sscanf(buf, fmt, &model.num))
			goto unknown;
		break;

	case MTRIX:
		if (0 > sscanf(buf, fmt, &mtrix.rowNum, &mtrix.serialNum,
				&mtrix.m1, &mtrix.m2, &mtrix.m3, &mtrix.v,
				&mtrix.given))
			goto unknown;
		break;

	case OBSLTE:
		if (0 > sscanf(buf, fmt, &obslte.continuation, obslte.timestamp,
				obslte.oldId, obslte.idMap[0],
				obslte.idMap[1], obslte.idMap[2],
				obslte.idMap[3], obslte.idMap[4],
				obslte.idMap[2], obslte.idMap[6],
				obslte.idMap[7]))
			goto unknown;
		break;

	case ORIGX:
		if (0 > sscanf(buf, fmt, &origx.rowNum, &origx.o1, &origx.o2,
				&origx.o3, &origx.t))
			goto unknown;
		break;

	case REVDAT:
		if (0 > sscanf(buf, fmt, &revdat.modification,
				&revdat.continuation, revdat.timestamp,
				revdat.id, &revdat.modType,
				revdat.corrections))
			goto unknown;
		break;

	case SCALE:
		if (0 > sscanf(buf, fmt, &scale.rowNum, &scale.s1, &scale.s2,
				&scale.s3, &scale.u))
			goto unknown;
		break;

	case SEQRES:
		if (0 > sscanf(buf, fmt, &seqres.serialNum, &seqres.chainId,
				&seqres.count, seqres.names[0], seqres.names[1],
				seqres.names[2], seqres.names[3],
				seqres.names[4], seqres.names[5],
				seqres.names[6], seqres.names[7],
				seqres.names[8], seqres.names[9],
				seqres.names[10], seqres.names[11],
				seqres.names[12]))
			goto unknown;
		break;

	case SHEET:
		sh = &sheet;
		sha0 = &sh->atoms[0].residue;
		sha1 = &sh->atoms[1].residue;
		if (0 > sscanf(buf, fmt, &sh->strandNum, sh->id, &sh->count,
				sh->residues[0].name, &sh->residues[0].chainId,
				&sh->residues[0].seqNum,
				&sh->residues[0].insertCode,
				sh->residues[1].name, &sh->residues[1].chainId,
				&sh->residues[1].seqNum,
				&sh->residues[1].insertCode, &sh->sense,
				sh->atoms[0].name, sha0->name, &sha0->chainId,
				&sha0->seqNum, &sha0->insertCode,
				sh->atoms[1].name, sha1->name, &sha1->chainId,
				&sha1->seqNum, &sha1->insertCode))
			goto unknown;
		break;

	case SITE:
		if (0 > sscanf(buf, fmt, &site.seqNum, site.id, &site.count,
				site.residues[0].name,
				&site.residues[0].chainId,
				&site.residues[0].seqNum,
				&site.residues[0].insertCode,
				site.residues[1].name,
				&site.residues[1].chainId,
				&site.residues[1].seqNum,
				&site.residues[1].insertCode,
				site.residues[2].name,
				&site.residues[2].chainId,
				&site.residues[2].seqNum,
				&site.residues[2].insertCode,
				site.residues[3].name,
				&site.residues[3].chainId,
				&site.residues[3].seqNum,
				&site.residues[3].insertCode))
			goto unknown;
		break;

	case SPRSDE:
		if (0 > sscanf(buf, fmt, &sprsde.continuation,
				sprsde.timestamp, sprsde.id,
				sprsde.supersede[0], sprsde.supersede[1],
				sprsde.supersede[2], sprsde.supersede[3],
				sprsde.supersede[4], sprsde.supersede[5],
				sprsde.supersede[6], sprsde.supersede[7]))
			goto unknown;
		break;

	case SSBOND:
		if (0 > sscanf(buf, fmt, &ssbond.seqNum,
				ssbond.residues[0].name,
				&ssbond.residues[0].chainId,
				&ssbond.residues[0].seqNum,
				&ssbond.residues[0].insertCode,
				ssbond.residues[1].name,
				&ssbond.residues[1].chainId,
				&ssbond.residues[1].seqNum,
				&ssbond.residues[1].insertCode,
				ssbond.comment))
			goto unknown;
		break;

	case SYMOP:
		if (0 > sscanf(buf, fmt, &symop.rowNum, &symop.serialNum,
				&symop.s1, &symop.s2, &symop.s3, &symop.t))
			goto unknown;
		break;

	case TER:
		if (0 > sscanf(buf, fmt, &ter.serialNum, ter.residue.name,
				&ter.residue.chainId, &ter.residue.seqNum,
				&ter.residue.insertCode))
			goto unknown;
		break;

	case TRNSFM:
		if (0 > sscanf(buf, fmt, &trnsfm.resultSerialNum,
				&trnsfm.applySerialNum,
				&trnsfm.sourceSerialNum))
			goto unknown;
		break;

	case TURN:
		if (0 > sscanf(buf, fmt, &turn.seqNum, turn.id,
				turn.residues[0].name,
				&turn.residues[0].chainId,
				&turn.residues[0].seqNum,
				&turn.residues[0].insertCode,
				turn.residues[1].name,
				&turn.residues[1].chainId,
				&turn.residues[1].seqNum,
				&turn.residues[1].insertCode, turn.comment))
			goto unknown;
		break;

	case TVECT:
		if (0 > sscanf(buf, fmt, &tvect.serialNum, &tvect.t1,
				&tvect.t2, &tvect.t3, tvect.comment))
			goto unknown;
		break;

user:
		rType = USER;
		fmt = pdbRecordFormatR[rType];
	case USER:
		if (0 > sscanf(buf, fmt, user.subtype, user.text))
			goto unknown;
		break;

	case USER_PDBRUN:
		if (0 > sscanf(buf, fmt, &userPdbrun.version))
			goto user;
		pdbrunInputVersion = userPdbrun.version;
		break;

	case USER_EYEPOS:
		if (0 > sscanf(buf, fmt, &userEyePos.xyz[0],
				&userEyePos.xyz[1], &userEyePos.xyz[2]))
			goto user;
		break;

	case USER_ATPOS:
		if (0 > sscanf(buf, fmt, &userAtPos.xyz[0],
				&userAtPos.xyz[1], &userAtPos.xyz[2]))
			goto user;
		break;

	case USER_WINDOW:
		if (0 > sscanf(buf, fmt, &userWindow.left, &userWindow.right,
				&userWindow.bottom, &userWindow.top,
				&userWindow.hither, &userWindow.yon))
			goto user;
		break;

	case USER_FOCUS:
		if (0 > sscanf(buf, fmt, &userFocus.focus))
			goto user;
		break;

	case USER_VIEWPORT:
		if (0 > sscanf(buf, fmt, &userViewport.xmin, &userViewport.xmax,
				&userViewport.ymin, &userViewport.ymax))
			goto user;
		break;

	case USER_BGCOLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userBgColor.rgb[0],
					&userBgColor.rgb[1],
					&userBgColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userBgColor.rgb[0],
				&userBgColor.rgb[1], &userBgColor.rgb[2]))
			goto user;
		break;

	case USER_ANGLE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userAngle.which,
					&userAngle.atom0, &userAngle.atom1,
					&userAngle.atom2, &userAngle.atom3,
					&userAngle.angle))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userAngle.atom0,
				&userAngle.atom1, &userAngle.atom2,
				&userAngle.atom3, &userAngle.angle))
			goto user;
		break;

	case USER_DISTANCE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userDistance.which,
					&userDistance.atom0,
					&userDistance.atom1,
					&userDistance.distance))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userDistance.atom0,
				&userDistance.atom1, &userDistance.distance))
			goto user;
		break;

	case USER_FILE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, userFile.filename))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userFile.model,
							userFile.filename))
			goto user;
		break;

	case USER_MARKNAME:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userMarkname.markname))
			goto user;
		break;

	case USER_MARK:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userMark.markname))
			goto user;
		break;

	case USER_CNAME:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userCName.name,
					&userCName.rgb[0], &userCName.rgb[1],
					&userCName.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userCName.rgb[0],
				&userCName.rgb[1], &userCName.rgb[2],
				userCName.name))
			goto user;
		break;

	case USER_COLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userColor.spec,
					&userColor.rgb[0], &userColor.rgb[1],
					&userColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userColor.rgb[0],
				&userColor.rgb[1], &userColor.rgb[2],
				userColor.spec))
			goto user;
		break;

	case USER_RADIUS:
		if (0 > sscanf(buf, fmt, &userRadius.radius))
			goto user;
		break;

	case USER_OBJECT:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userObject.model))
				goto user;
		}
		break;

	case USER_ENDOBJ:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userEndObj.model))
				goto user;
		}
		break;

	case USER_CHAIN:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userChain.atom0,
					&userChain.atom1))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userChain.atom0,
				&userChain.atom1))
			goto user;
		break;

	case USER_GFX_BEGIN:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userGfxBegin.unknown))
			goto user;
		userGfxBegin.primitive = getGfxType(userGfxBegin.unknown);
		break;

	case USER_GFX_END:
		if (pdbrunInputVersion < 6)
			goto user;
		break;

	case USER_GFX_COLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userGfxColor.spec,
					&userGfxColor.rgb[0],
					&userGfxColor.rgb[1],
					&userGfxColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userGfxColor.rgb[0],
				&userGfxColor.rgb[1], &userGfxColor.rgb[2],
				userGfxColor.spec))
			goto user;
		break;

	case USER_GFX_NORMAL:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxNormal.xyz[0],
				&userGfxNormal.xyz[1],
				&userGfxNormal.xyz[2]))
			goto user;
		break;

	case USER_GFX_VERTEX:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxVertex.xyz[0],
				&userGfxVertex.xyz[1],
				&userGfxVertex.xyz[2]))
			goto user;
		break;

	case USER_GFX_FONT:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userGfxFont.name,
					&userGfxFont.size))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userGfxFont.size,
				userGfxFont.name))
			goto user;
		break;

	case USER_GFX_TEXTPOS:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxTextPos.xyz[0],
				&userGfxTextPos.xyz[1], &userGfxTextPos.xyz[2]))
			goto user;
		break;

	case USER_GFX_LABEL:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userGfxLabel.xyz[0],
					&userGfxLabel.xyz[1],
					&userGfxLabel.xyz[2],
					userGfxLabel.text))
				goto user;
		} else if (0 > sscanf(buf, fmt, userGfxLabel.text))
			goto user;
		// TODO: process text?
		break;

	case USER_GFX_MOVE:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxMove.xyz[0],
				&userGfxMove.xyz[1], &userGfxMove.xyz[2]))
			goto user;
		break;

	case USER_GFX_DRAW:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxDraw.xyz[0],
				&userGfxDraw.xyz[1], &userGfxDraw.xyz[2]))
			goto user;
		break;

	case USER_GFX_MARKER:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxMarker.xyz[0],
				&userGfxMarker.xyz[1],
				&userGfxMarker.xyz[2]))
			goto user;
		break;

	case USER_GFX_POINT:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxPoint.xyz[0],
				&userGfxPoint.xyz[1], &userGfxPoint.xyz[2]))
			goto user;
		break;
	}
}

} // namespace PDBPP

//
//	Copyright (c) 1989,1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_sprntf.cc,v 1.4 94/09/06 15:03:07 gregc Exp $
//


#include <ctype.h>
#include <cstring>
#include <cstdarg>

#define	OVERFLOW_CHAR	'*'

namespace PDBPP {

using namespace std;

	// scratch must be big enough to hold the largest number
static char	scratch[1024];


static char	*outint(int, int, int, char, char, int, char *, char);
static char	*outunsigned(unsigned int, int, char, int, char *);
static char	*outstr(char *, int, int, char, int, char *);
static char	*outfloat(double, int, int, char, int, char *);
static char	*outexp(double, int, int, char, int, char *);
static char	*e_out(int, char *);

int
PDB::sprintf(char *outbuf, const char *fmt, ...)
{
	va_list		argv;
	char		*p;
	const char	*f;
	int		field1, field2;
	char		c, fill_char;
	int		inum;
	unsigned 	unum;
	double		fnum;
	int		left_justify;

	va_start(argv, fmt);
	f = fmt;
	p = outbuf;
	while (*f) {
		if (*f == '%') {
			f++;
			if (*f == '-')
				left_justify = 1, f++;
			else
				left_justify = 0;

			if (*f == '0')
				fill_char = '0', f++;
			else
				fill_char = ' ';

			if (isdigit(*f)) {
				field1 = *f++ - '0';
				while (isdigit(*f))
					field1 = field1 * 10 + *f++ - '0';
			}
			else
				field1 = -1;

			if (*f == '.') {
				f++;
				field2 = 0;
				while (isdigit(*f))
					field2 = field2 * 10 + *f++ - '0';
			}
			else
				field2 = -1;

			if (*f == 'l' || *f == 'h')
				f++;

			while (isspace(*f))
				f++;
			switch (*f) {
			  case 'c':
				c = (char) va_arg(argv, int);
				if (c == '\0')
					c = ' ';
				if (left_justify)
					*p++ = c;
				while (--field1 > 0)
					*p++ = fill_char;
				if (!left_justify)
					*p++ = c;
				break;
			  case 'd':
			  case 'D':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 10, fill_char, 'a',
					left_justify, p, (*f == 'D') ? ' ':'0');
				break;
			  case 'e':
				fnum = va_arg(argv, double);
				if (field2 < 0)
					field2 = 6;
				p = outexp(fnum, field1, field2, fill_char,
					left_justify, p);
				break;
			  case 'f':
				fnum = va_arg(argv, double);
				if (field2 < 0)
					field2 = 6;
				p = outfloat(fnum, field1, field2, fill_char,
					left_justify, p);
				break;
			  case 'o':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 8, fill_char, 'a',
					left_justify, p, '0');
				break;
			  case 's':
				p = outstr(va_arg(argv, char *), field1, field2,
					fill_char, left_justify, p);
				break;
			  case 'u':
				unum = va_arg(argv, unsigned);
				p = outunsigned(unum, field1, fill_char,
					left_justify, p);
				break;
			  case 'x':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 16, fill_char, 'a',
					left_justify, p, '0');
				break;
			  case 'X':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 16, fill_char, 'A',
					left_justify, p, '0');
				break;
			  default:
				if (left_justify)
					*p++ = *f;
				while (--field1 > 0)
					*p++ = fill_char;
				if (!left_justify)
					*p++ = *f;
				break;
			}
			f++;
		}
		else if (*f == '\\') {		/* Special character */
			switch (*++f) {
			  case 'n':
				*p++ = '\n';
				break;
			  case 'r':
				*p++ = '\r';
				break;
			  case 'b':
				*p++ = '\b';
				break;
			  case 't':
				*p++ = '\t';
				break;
			  case 'f':
				*p++ = '\f';
				break;
			  case '0': case '1': case '2': case '3':
			  case '4': case '5': case '6': case '7':
				inum = *f++ - '0';
				if (*f >= '0' && *f <= '7') {
					inum = inum * 8 + *f++ - '0';
					if (*f >= '0' && *f <= '7')
						inum = inum * 8 + *f++ - '0';
				}
				f--;
				*p++ = (char) inum;
				break;
			  default:
				*p++ = *f;
			}
			f++;
		}
		else				/* Normal character */
			*p++ = *f++;
	}
	*p = '\0';
	va_end(argv);
	return p - outbuf;
}

static char *
e_out(int width, char *where)
{
	while (width-- > 0)
		*where++ = OVERFLOW_CHAR;
	return where;
}

static char *
outint(int value, int width, int radix, char fill_char, char hex,
					int left_justify, char *p, char zero)
{
	char	*s;
	int	n;
	int	negative;

	if (value < 0)
		negative = 1, value = -value, width--;
	else
		negative = 0;
	s = scratch;
	if (value)
		do {
			n = value % radix;
			*s++ = n < 10 ? '0' + n : hex + n - 10;
			value /= radix;
		} while (value);
	else
		*s++ = zero;
	n = s - scratch;
	if (width != -1 && n > width)
		return e_out(width + negative, p);

	if (negative && fill_char == '0')
		*p++ = '-';
	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	if (negative && fill_char == ' ')
		*p++ = '-';
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}

static char *
outunsigned(unsigned int value, int width, char fill_char, int left_justify,
									char *p)
{
	char	*s;
	int	n;

	s = scratch;
	while (value) {
		*s++ = value % 10 + '0';
		value /= 10;
	}
	n = s - scratch;
	if (n == 0)
		*s++ = '0', n = 1;
	if (width != -1 && n > width)
		return e_out(width, p);

	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}

static char *
outstr(char *s, int width, int maxstr, char fill_char, int left_justify, char *p)
{
	int	len;

	len = strlen(s);
	if (maxstr >= 0 && len > maxstr)
		len = maxstr;
	if (width != -1 && len > width)
		return e_out(width, p);

	if (!left_justify)
		while (width-- > len)
			*p++ = fill_char;
	else
		width -= len;
	while (len--)
		*p++ = *s++;
	if (left_justify)
		while (width-- > 0)
			*p++ = fill_char;
	return p;
}

static char *
outfloat(double value, int width, int nplace, char fill_char, int left_justify,
									char *p)
{
	int	i, intval;
	char	*place, *to, *from;
	int	negative;

	negative = value < 0.0 ? 1 : 0;
		
	if (negative)
		value = -value;

	for (i = 0; i < nplace; i++)
		value *= 10.0;

	intval = (int) (value + 0.5);

	if (width == -1)
		width = nplace + 4;		/* TODO: fix */
	else if (nplace + (nplace == 0 ? 1 : 2) > width)
		return e_out(width, p);

	for (place = p + width - 1; place >= p + width - nplace; place--) {
		*place = '0' + intval % 10;
		intval /= 10;
	}

	if (nplace > 0)
		*place-- = '.';

	if (intval == 0)
		*place-- = '0';

	for (; place >= p; place--) {
		if (intval == 0)
			break;
		else {
			*place = '0' + intval % 10;
			intval /= 10;
		}
	}

	if (intval != 0)
		return e_out(width, p);

	if (place < p && negative)
		return e_out(width, p);

	if (left_justify) {
		for (from = place + 1, to = (negative ? p + 1 : p);
						from < p + width; from++, to++)
			*to = *from;
		for (; to < p + width; to++)
			*to = fill_char;
		if (negative)
			*p = '-';
	} else {
		for (to = place; to >= p; to--)
			*to = fill_char;
		if (negative)
			if (fill_char == ' ')
				*place = '-';
			else
				*p = '-';
	}

	return p + width;
}

static char *
outexp(double value, int width, int nplace, char fill_char, int left_justify,
									char *p)
{
	int	n;
	char	*s;
	int	negative;
	double	fraction;

	if (value < 0)
		negative = 1, value = -value, width--;
	else
		negative = 0;

	n = 0;
	while (value > 10)
		n++, value /= 10;
	if (value)
		while (value < 1)
			n--, value *= 10;

	s = scratch;
	if (n < 0) {
		n = -n;
		*s++ = n % 10 + '0';
		*s++ = n / 10 + '0';
		*s++ = '-';
	}
	else {
		*s++ = n % 10 + '0';
		*s++ = n / 10 + '0';
		*s++ = '+';
	}
	*s = 'e';

	s = scratch + nplace + 4;	/* 4 == strlen("e+00") */
	fraction = value - (int) value;
	for (n = 0; n < nplace; n++) {
		fraction *= 10.0;
		*--s = '0' + (int) fraction;
		fraction -= (int) fraction;
	}

	s = scratch + nplace + 4;
	if (nplace)
		*s++ = '.';
	n = (int) value;
	if (n)
		*s++ = n % 10 + '0';
	else
		*s++ = '0';
	n = s - scratch;
	if (width != -1 && n > width)
		return e_out(width + negative, p);

	if (negative && fill_char == '0')
		*p++ = '-';
	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	if (negative && fill_char == ' ')
		*p++ = '-';
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}

} // namespace PDBPP
//
//	Copyright (c) 1989,1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_sscanf.cc,v 1.8 94/09/06 15:03:23 gregc Exp $
//


#include <ctype.h>
#include <cstdarg>
#include <cstdlib>

//
//	pdb_sscanf performs similarly to sscanf, execept that fields are of
//	fixed length and a complete line is always consumed.  The field
//	width defaults to one.  If the line is shorter than expected then
//	the default is returned.
//
//		d	get an integer.  Default:  0.
//		f	get a floating point number (C double).  Default:  0.0.
//		(space) ignore characters within field
//		s	get a C string, leading and trailing spaces are
//			stripped; the field width is used as a limit on
//			the string length, the null character is appended
//			to the end of the string.  Default:  empty string.
//		c	get a character(s); no stripping of spaces, nor is
//			a null character appended.  Default:  space(s).
//

#define	MAXFIELDSIZE	264

namespace PDBPP {

int
PDB::sscanf(const char *buffer, const char *fmt, ...)
{
    //DEBUG:
    //std::cout << "PDB::sscanf: " << buffer << " : " << strlen(buffer) << " :  " << fmt << "\n";
	va_list	ap;
	int	i, field_width;
	int	nmatch;
	char	*s, *t;
	char	tmp[MAXFIELDSIZE];

	va_start(ap, fmt);
	nmatch = 0;
	for (; *fmt != '\0'; fmt++) {
		if (*fmt != '%') {
			if (*buffer == *fmt)
				buffer++;
			else if (*buffer != '\0' && *buffer != '\n')
				return -1;
			continue;
		}

		// calculate field_width
		field_width = 0;
		for (++fmt; isdigit(*fmt); fmt++)
			field_width = field_width * 10 + *fmt - '0';
		if (field_width == 0)
			field_width = 1;	// default
		if (*buffer != '\0' && *buffer != '\n')
			nmatch++;

		switch (*fmt) {

		case 'd':			// integer
			// if we've already seen the end of the buffer, don't
			// try to get anymore characters
			if (*buffer == '\0' || *buffer == '\n') {
				*(va_arg(ap, int *)) = 0;
				break;
			}

			s = tmp;
			for (i = 0; i < field_width; i++) {
				if (*buffer == '\0' || *buffer == '\n')
					break;
				*s++ = *buffer++;
			}
			*s = '\0';
			// remove trailing spaces
			while (s > tmp && isspace(*(s - 1)))
				*--s = '\0';
			*(va_arg(ap, int *)) = (int) strtol(tmp, &t, 10);
			if (t != s)
				return -1;
			break;

		case 'f':			// floating point
			// if we've already seen the end of the buffer, don't
			// try to get anymore characters
			if (*buffer == '\0' || *buffer == '\n') {
				*(va_arg(ap, double *)) = 0.0;
				break;
			}

			s = tmp;
			for (i = 0; i < field_width; i++) {
				if (*buffer == '\0' || *buffer == '\n')
					break;
				*s++ = *buffer++;
			}
			*s = '\0';
			// remove trailing spaces
			while (s > tmp && isspace(*(s - 1)))
				*--s = '\0';
			*(va_arg(ap, double *)) = strtod(tmp, &t);
			if (t != s)
				return -1;
			break;

		case 's':			// string
			// if we've already seen the end of the buffer, don't
			// try to get anymore characters
			if (*buffer == '\0' || *buffer == '\n') {
				*(va_arg(ap, char *)) = '\0';
				break;
			}

			s = t = va_arg(ap, char *);
			for (i = 0; i < field_width; i++) {
				if (*buffer == '\0' || *buffer == '\n')
					break;
				*s++ = *buffer++;
			}
			*s = '\0';
			// remove trailing spaces
			while (s > t && isspace(*--s))
				*s = '\0';
			break;

		case 'c':			// character(s)
			s = va_arg(ap, char *);
			for (i = 0; i < field_width; i++)
				s[i] = ' ';	// default

			// if we've already seen the end of the buffer, don't
			// try to get anymore characters
			if (*buffer == '\0' || *buffer == '\n')
				break;

			for (i = 0; i < field_width; i++) {
				if (*buffer == '\0' || *buffer == '\n')
					break;
				*s++ = *buffer++;
			}
			break;

		case ' ':			// space (ignore)
			// if we've already seen the end of the buffer, don't
			// try to get anymore characters
			if (*buffer == '\0' || *buffer == '\n')
				break;

			for (i = 0; i < field_width; i++, buffer++)
				if (*buffer == '\0' || *buffer == '\n')
					break;
			break;

		default:
			va_end(ap);
			return -1;
		}
	}
	va_end(ap);
	return nmatch;
}

} // namespace PDBPP
//
//	Copyright (c) 1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_type.cc,v 1.14 1995/05/22 19:53:09 gregc Exp $
//
//	subroutine for reading PDB format files
//


# include	<ctype.h>
# include	<cstring>

#if 0
extern "C" int strcasecmp (const char *s1, const char *s2);
extern "C" int strncasecmp (const char *s1, const char *s2, size_t n);
#endif

# ifndef _toupper
# define	_toupper	toupper
# endif

namespace PDBPP {

using namespace std;

int PDB::pdbrunInputVersion = PDB::PDBRUNVersion;	// just in case
int PDB::pdbrunOutputVersion = PDB::PDBRUNVersion;	// just in case

PDB::GfxType
PDB::getGfxType(const char *buf)
{
	switch (buf[0]) {
	case 'L': case 'l':
		if (strcasecmp(buf + 1, "INE-LOOP") == 0)
			return GFX_LINE_LOOP;
		if (strcasecmp(buf + 1, "INE-STRIP") == 0)
			return GFX_LINE_STRIP;
		if (strcasecmp(buf + 1, "INES") == 0)
			return GFX_LINES;
		break;
	case 'M': case 'm':
		if (strcasecmp(buf + 1, "ARKERS") == 0)
			return GFX_MARKERS;
		break;
	case 'P': case 'p':
		if (strcasecmp(buf + 1, "OINTS") == 0)
			return GFX_POINTS;
		if (strcasecmp(buf + 1, "OLYGON") == 0)
			return GFX_POLYGON;
		break;
	case 'Q': case 'q':
		if (strcasecmp(buf + 1, "UAD-STRIP") == 0)
			return GFX_QUAD_STRIP;
		if (strcasecmp(buf + 1, "UADS") == 0)
			return GFX_QUADS;
		break;
	case 'T': case 't':
		if (strcasecmp(buf + 1, "RIANGLE-FAN") == 0)
			return GFX_TRIANGLE_FAN;
		if (strcasecmp(buf + 1, "RIANGLE-STRIP") == 0)
			return GFX_TRIANGLE_STRIP;
		if (strcasecmp(buf + 1, "RIANGLES") == 0)
			return GFX_TRIANGLES;
		break;
	}
	return GFX_UNKNOWN;
}

const char *
PDB::gfxChars(GfxType gt)
{
	switch (gt) {
	default:		return "UNKNOWN";
	case GFX_POINTS:	return "POINTS";
	case GFX_MARKERS:	return "MARKERS";
	case GFX_LINES:		return "LINES";
	case GFX_LINE_STRIP:	return "LINE-STRIP";
	case GFX_LINE_LOOP:	return "LINE-LOOP";
	case GFX_TRIANGLES:	return "TRIANGLES";
	case GFX_TRIANGLE_STRIP:	return "TRIANGLE-STRIP";
	case GFX_TRIANGLE_FAN:	return "TRIANGLE-FAN";
	case GFX_QUADS:		return "QUADS";
	case GFX_QUAD_STRIP:	return "QUAD-STRIP";
	case GFX_POLYGON:	return "POLYGON";
	}
}

static PDB::RecordType
pdbrun5Type(const char *buf)
{
	switch (buf[0]) {
	case 'A': case 'a':
		if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
			return PDB::USER_ANGLE;
		if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
			return PDB::USER_ATPOS;
		break;
	case 'B': case 'b':
		if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
			return PDB::USER_BGCOLOR;
		break;
	case 'C': case 'c':
		if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
			return PDB::USER_CHAIN;
		if (strncasecmp(buf + 1, "NAME ", 5) == 0)
			return PDB::USER_CNAME;
		if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
			return PDB::USER_COLOR;
		break;
	case 'D': case 'd':
		if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
			return PDB::USER_DISTANCE;
		break;
	case 'E': case 'e':
		if (strncasecmp(buf + 1, "NDOBJ ", 6) == 0)
			return PDB::USER_ENDOBJ;
		if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
			return PDB::USER_EYEPOS;
		break;
	case 'F': case 'f':
		if (strncasecmp(buf + 1, "ILE ", 4) == 0)
			return PDB::USER_FILE;
		if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
			return PDB::USER_FOCUS;
		break;
	case 'G': case 'g':
		if (strncasecmp(buf + 1, "FX ", 3) != 0)
			break;
		if (strncasecmp(buf + 4, "COLOR ", 6) == 0)
			return PDB::USER_GFX_COLOR;
		if (strncasecmp(buf + 4, "DRAW ", 5) == 0)
			return PDB::USER_GFX_DRAW;
		if (strncasecmp(buf + 4, "FONT ", 5) == 0)
			return PDB::USER_GFX_FONT;
		if (strncasecmp(buf + 4, "LABEL ", 6) == 0)
			return PDB::USER_GFX_LABEL;
		if (strncasecmp(buf + 4, "MARKER ", 7) == 0)
			return PDB::USER_GFX_MARKER;
		if (strncasecmp(buf + 4, "MOVE ", 5) == 0)
			return PDB::USER_GFX_MOVE;
		if (strncasecmp(buf + 4, "POINT ", 6) == 0)
			return PDB::USER_GFX_POINT;
		break;
	case 'O': case 'o':
		if (strncasecmp(buf + 1, "BJECT ", 6) == 0)
			return PDB::USER_OBJECT;
		break;
	case 'P': case 'p':
		if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
			return PDB::USER_PDBRUN;
		break;
	case 'R': case 'r':
		if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
			return PDB::USER_RADIUS;
		break;
	case 'V': case 'v':
		if (strncasecmp(buf + 1, "IEWPORT ", 8) == 0)
			return PDB::USER_VIEWPORT;
		break;
	case 'W': case 'w':
		if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
			return PDB::USER_WINDOW;
		break;
	}
	return PDB::USER;
}

static PDB::RecordType
pdbrun6Type(const char *buf)
{
	switch (buf[0]) {
	case 'A': case 'a':
		if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
			return PDB::USER_ANGLE;
		if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
			return PDB::USER_ATPOS;
		break;
	case 'B': case 'b':
		if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
			return PDB::USER_BGCOLOR;
		break;
	case 'C': case 'c':
		if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
			return PDB::USER_CHAIN;
		if (strncasecmp(buf + 1, "NAME ", 5) == 0)
			return PDB::USER_CNAME;
		if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
			return PDB::USER_COLOR;
		break;
	case 'D': case 'd':
		if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
			return PDB::USER_DISTANCE;
		break;
	case 'E': case 'e':
		if (strncasecmp(buf + 1, "NDOBJ", 5) == 0
		&& (buf[6] == '\0' || buf[6] == '\n' || buf[6] == ' '))
			return PDB::USER_ENDOBJ;
		if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
			return PDB::USER_EYEPOS;
		break;
	case 'F': case 'f':
		if (strncasecmp(buf + 1, "ILE ", 4) == 0)
			return PDB::USER_FILE;
		if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
			return PDB::USER_FOCUS;
		break;
	case 'G': case 'g':
		if (buf[1] != 'F' || buf[2] != 'X' || buf[3] != ' ')
			break;
		switch (buf[4]) {
		case 'B': case 'b':
			if (strncasecmp(buf + 5, "EGIN ", 5) == 0)
				return PDB::USER_GFX_BEGIN;
			break;
		case 'C': case 'c':
			if (strncasecmp(buf + 5, "OLOR ", 5) == 0)
				return PDB::USER_GFX_COLOR;
			break;
		case 'E': case 'e':
			if (buf[5] == 'N' && buf[6] == 'D'
			&& (buf[7] == '\0' || buf[7] == '\n' || buf[7] == ' '))
				return PDB::USER_GFX_END;
			break;
		case 'F': case 'f':
			if (strncasecmp(buf + 5, "ONT ", 4) == 0)
				return PDB::USER_GFX_FONT;
			break;
		case 'L': case 'l':
			if (strncasecmp(buf + 5, "ABEL ", 5) == 0)
				return PDB::USER_GFX_LABEL;
			break;
		case 'N': case 'n':
			if (strncasecmp(buf + 5, "ORMAL ", 6) == 0)
				return PDB::USER_GFX_NORMAL;
			break;
		case 'T': case 't':
			if (strncasecmp(buf + 5, "EXTPOS ", 7) == 0)
				return PDB::USER_GFX_TEXTPOS;
			break;
		case 'V': case 'v':
			if (strncasecmp(buf + 5, "ERTEX ", 6) == 0)
				return PDB::USER_GFX_VERTEX;
			break;
		}
		break;
	case 'M': case 'm':
		if (strncasecmp(buf + 1, "ARK ", 4) == 0)
			return PDB::USER_MARK;
		if (strncasecmp(buf + 1, "ARKNAME ", 6) == 0)
			return PDB::USER_MARKNAME;
		break;
	case 'O': case 'o':
		if (strncasecmp(buf + 1, "BJECT", 5) == 0
		&& (buf[6] == '\0' || buf[6] == '\n' || buf[6] == ' '))
			return PDB::USER_OBJECT;
		break;
	case 'P': case 'p':
		if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
			return PDB::USER_PDBRUN;
		break;
	case 'R': case 'r':
		if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
			return PDB::USER_RADIUS;
		break;
	case 'V': case 'v':
		if (strncasecmp(buf + 1, "IEWPORT ", 6) == 0)
			return PDB::USER_VIEWPORT;
		break;
	case 'W': case 'w':
		if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
			return PDB::USER_WINDOW;
		break;
	}
	return PDB::USER;
}

PDB::RecordType
PDB::getType(const char *buf)
{
    //DEBUG:
    //std::cout << buf << '\n';

	char	rt[4];		// PDB record type
	int	i;

	for (i = 0; buf[i] != '\0' && buf[i] != '\n' && i < 4; i += 1) {
		if (islower(buf[i]))
			rt[i] = _toupper(buf[i]);
		else
			rt[i] = buf[i];
	}
	
    //DEBUG:
    //std::cout << rt << '\n';
	if (i < 4)
		for (; i < 4; i += 1)
			rt[i] = ' ';

	switch (rt[0]) {

	case 'A':
		switch (rt[1]) {
		case 'G':
			if (rt[2] == 'R' && rt[3] == 'D')
				return AGRDES;
			if (rt[2] == 'G' && rt[3] == 'R')
				return AGGRGT;
			break;
		case 'N':
			if (rt[2] == 'I' && rt[3] == 'S')
				return ANISOU;
			break;
		case 'T':
			if (rt[2] == 'O' && rt[3] == 'M') {
			    //DEBUG:
			    //std::cout << "ATOM\n";
				return ATOM;
				}
			break;
		case 'U':
			if (rt[2] == 'T' && rt[3] == 'H')
				return AUTHOR;
			break;
		}
		break;

	case 'C':
		switch (rt[1]) {
		case 'M':
			if (rt[2] == 'P' && rt[3] == 'D')
				return CMPDES;
			if (rt[2] == 'P' && rt[3] == 'O')
				return CMPONT;
			break;
		case 'O':
			if (rt[2] == 'M' && rt[3] == 'P')
				return COMPND;
			if (rt[2] == 'N' && rt[3] == 'E')
				return CONECT;
			break;
		case 'R':
			if (rt[2] == 'Y' && rt[3] == 'S')
				return CRYST1;
			break;
		}
		break;

	case 'E':
		switch (rt[1]) {
		case 'N':
			if (rt[2] == 'D' && rt[3] == ' ')
				return END;
			if (rt[2] == 'D' && rt[3] == 'M')
				return ENDMDL;
			break;
		case 'X':
			if (rt[2] == 'P' && rt[3] == 'D')
				return EXPDTA;
			break;
		}
		break;

	case 'F':
		switch (rt[1]) {
		case 'T':
			if (rt[2] == 'N' && rt[3] == 'O')
				return FTNOTE;
			break;
		case 'O':
			if (rt[2] == 'R' && rt[3] == 'M')
				return FORMUL;
			break;
		}
		break;

	case 'H':
		if (rt[1] != 'E')
			break;
		if (rt[2] == 'T' && rt[3] == 'A')
			return HETATM;
		if (rt[2] == 'A' && rt[3] == 'D')
			return HEADER;
		if (rt[2] == 'T' && rt[3] == ' ')
			return HET;
		if (rt[2] == 'L' && rt[3] == 'I')
			return HELIX;
		break;

	case 'J':
		if (rt[1] == 'R' && rt[2] == 'N' && rt[3] == 'L')
			return JRNL;
		break;

	case 'M':
		switch (rt[1]) {
		case 'A':
			if (rt[2] == 'S' && rt[3] == 'T')
				return MASTER;
			break;
		case 'O':
			if (rt[2] == 'D' && rt[3] == 'E')
				return MODEL;
			break;
		case 'T':
			if (rt[2] == 'R' && rt[3] == 'I')
				return MTRIX;
			if (rt[2] == 'X' && rt[3] == 'D')
				return MTXDES;
			break;
		}
		break;

	case 'O':
		switch (rt[1]) {
		case 'B':
			if (rt[2] == 'S' && rt[3] == 'L')
				return OBSLTE;
			break;
		case 'R':
			if (rt[2] == 'I' && rt[3] == 'G')
				return ORIGX;
			break;
		}
		break;

	case 'R':
		if (rt[1] != 'E')
			break;
		if (rt[2] == 'M' && rt[3] == 'A')
			return REMARK;
		if (rt[2] == 'V' && rt[3] == 'D')
			return REVDAT;
		break;

	case 'S':
		switch (rt[1]) {

		case 'C':
			if (rt[2] == 'A' && rt[3] == 'L')
				return SCALE;
			break;

		case 'E':
			if (rt[2] == 'Q' && rt[3] == 'R')
				return SEQRES;
			break;

		case 'H':
			if (rt[2] == 'E' && rt[3] == 'E')
				return SHEET;
			break;

		case 'I':
			if (rt[2] == 'T' && rt[3] == 'E')
				return SITE;
			if (rt[2] == 'G' && rt[3] == 'A')
				return SIGATM;
			if (rt[2] == 'G' && rt[3] == 'U')
				return SIGUIJ;
			break;

		case 'O':
			if (rt[2] == 'U' && rt[3] == 'R')
				return SOURCE;
			break;

		case 'P':
			if (rt[2] == 'R' && rt[3] == 'S')
				return SPRSDE;
			break;

		case 'S':
			if (rt[2] == 'B' && rt[3] == 'O')
				return SSBOND;
			break;

		case 'Y':
			if (rt[2] == 'M' && rt[3] == 'D')
				return SYMDES;
			if (rt[2] == 'M' && rt[3] == 'O')
				return SYMOP;
			break;
		}
		break;

	case 'T':
		switch (rt[1]) {
		case 'E':
			if (rt[2] == 'R' && rt[3] == ' ')
				return TER;
			break;
		case 'R':
			if (rt[2] == 'N' && rt[3] == 'S')
				return TRNSFM;
			break;
		case 'U':
			if (rt[2] == 'R' && rt[3] == 'N')
				return TURN;
			break;
		case 'V':
			if (rt[2] == 'E' && rt[3] == 'C')
				return TVECT;
			break;
		}
		break;

	case 'U':
		if (rt[1] == 'S' && rt[2] == 'E' && rt[3] == 'R')
			switch (pdbrunInputVersion) {
			case 1: case 2: case 3: case 4: case 5:
				return pdbrun5Type(buf + 6);
			case 6:
				return pdbrun6Type(buf + 6);
			default:
				if (strncasecmp(buf + 6, "PDBRUN ", 7) == 0)
					return USER_PDBRUN;
				return USER;
			}
		break;
	}
	return UNKNOWN;
}

} // namespace PDBPP
