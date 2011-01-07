#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "trace_controller.hh"
#include "multiple_alignment.hh"
#include "sequence.hh"

#include <math.h>
#include <assert.h>

using namespace locarna;


TraceController::TraceRange::seqentry_pair_t
TraceController::TraceRange::
remove_common_gaps(const MultipleAlignment::SeqEntry &aliA,
		   const MultipleAlignment::SeqEntry &aliB) {
    size_t lenAli = aliA.seq().length();
    
    std::string raliA="";
    std::string raliB="";

    for (size_t i=1; i<=lenAli;i++) {
	if (!(MultipleAlignment::SeqEntry::is_gap_symbol(aliA.seq()[i])
	      && MultipleAlignment::SeqEntry::is_gap_symbol(aliB.seq()[i])
	      )) {
	    raliA+=aliA.seq()[i];
	    raliB+=aliB.seq()[i];
	}
    }
    
    return 
	seqentry_pair_t(
			MultipleAlignment::SeqEntry("raliA",raliA),
			MultipleAlignment::SeqEntry("raliB",raliB)
			);
}   


TraceController::TraceRange
::TraceRange(
	     const MultipleAlignment::SeqEntry &pseqA,
	     const MultipleAlignment::SeqEntry &pseqB,
	     const MultipleAlignment::SeqEntry &paliA,
	     const MultipleAlignment::SeqEntry &paliB,
	     size_type delta) {
    
    // pseqA and pseqB can contain gaps, therefore we call these strings profile sequences    

    // std::cout << pseqA.seq().to_string() << std::endl
    // 	      << paliA.seq().to_string() << std::endl
    // 	      << paliB.seq().to_string() << std::endl
    // 	      << pseqB.seq().to_string() << std::endl; 
	
    assert(paliA.seq().length() == paliB.seq().length());
    
    size_t plenA = pseqA.seq().length();
    size_t plenB = pseqB.seq().length();
    
    min_col_vector.resize(plenA+1);
    max_col_vector.resize(plenA+1);
    
    seqentry_pair_t ali =
	remove_common_gaps(paliA,paliB);
    
    const MultipleAlignment::SeqEntry &aliA=ali.first;
    const MultipleAlignment::SeqEntry &aliB=ali.second;
      
    size_t lenAli = aliA.seq().length();
    size_t lenA = pseqA.length_wogaps();
    size_t lenB = pseqB.length_wogaps();
        
    // std::cout << plenA << " "
    // 	      << plenB << " "
    // 	      << lenA << " "
    // 	      << lenB << " "
    // 	      << lenAli << " "
    // 	      << std::endl;
    
#ifdef COLUMN_CUT_DISTANCE
    // this code will compute the permissible cuts according to a definition
    // of the alignment deviation that limits the column cut distance to Delta
    
    // iterate over positions in sequence A
    for (size_t pi=0; pi <= plenA; pi++) {
	
	size_t left_i=pseqA.col_to_pos(pi).first; // left_i is the position in seqA left of the gap 
	size_t right_i=pseqA.col_to_pos(pi+1).second; // right_i is a position in seqA right of the gap
	// a gap starting at position pi in plenA corresponds to a gap between positions left_i, right_i in seqA
	
	size_t left_col = aliA.pos_to_col(left_i);
	size_t right_col = aliA.pos_to_col(right_i);
	
	// add delta deviation to columns
	left_col = std::max(delta,left_col)-delta;
	right_col = std::min(lenAli+1,right_col+delta);
	
	size_t left_j = aliB.col_to_pos(left_col).first;
	size_t right_j = aliB.col_to_pos(right_col).second;
	
	size_t left_pj = pseqB.pos_to_col(left_j);
	size_t right_pj = pseqB.pos_to_col(right_j);
	
	
	min_col_vector[pi] = left_pj;
	max_col_vector[pi] = right_pj-1;

    }
#else // POSITION_CUT_DISTANCE
    // this code computes the permissible cuts according to a definition
    // of the alignment deviation that limits the position cut distance to Delta
    
    // initialize col vectors
    for (size_t pi=0; pi <= plenA; pi++) {
	min_col_vector[pi] = plenB;
	max_col_vector[pi] = 0;
    }
    
    // iterate over columns of the alignment paliA/paliB
    for (size_t c=0; c <= lenAli; c++) {
	
	// determine cut ^t(pi,pj) of the alignment ^t(paliA,paliB) at column c
	size_t i = aliA.col_to_pos(c).first; // position in sequence A that corresponds to column c
	size_t j = aliB.col_to_pos(c).first; // position in sequence B that corresponds to column c
	
	// std::cout << c << " "
	// 	  << i << " "
	// 	  << j << " "
	// 	  << std::endl; 
	
	// this cut corresponds to a set of cuts C in alignments of pseqA and pseqB,
	// we describe this set by two ranges pi_min..pi_max and pj_min..pj_max
	size_t pi_min=pseqA.pos_to_col(i);
	size_t pi_max=pseqA.pos_to_col(i+1)-1;
	size_t pj_min=pseqB.pos_to_col(j);
	size_t pj_max=pseqB.pos_to_col(j+1)-1;

	// std::cout <<"  "
	// 	  << pi_min << " "
	// 	  << pi_max << " "
	// 	  << pj_min << " "
	// 	  << pj_max << " "
	// 	  << std::endl; 
	
	// determine the positions in delta distance 
	size_t i_minus = std::max(delta,i)-delta;
	size_t i_plus  = std::min(lenA,i+delta);
	size_t j_minus = std::max(delta,j)-delta;
	size_t j_plus  = std::min(lenB,j+delta);

	// std::cout <<"  "
	// 	  << i_minus << " "
	// 	  << i_plus << " "
	// 	  << j_minus << " "
	// 	  << j_plus << " "
	// 	  << std::endl; 

	
	// project to positions in pseqA and pseqB respectively
	size_t pi_minus = pseqA.pos_to_col(i_minus);
	size_t pi_plus  = pseqA.pos_to_col(i_plus);
	size_t pj_minus = pseqB.pos_to_col(j_minus);
	size_t pj_plus  = pseqB.pos_to_col(j_plus);
	
	// std::cout <<"  "
	// 	  << pi_minus << " "
	// 	  << pi_plus << " "
	// 	  << pj_minus << " "
	// 	  << pj_plus << " "
	// 	  << std::endl;


	// for all cuts in C, potentially update min_col_vector and max_col_vector
	for (size_t pi=pi_min; pi<=pi_max; pi++) {
	    min_col_vector[pi] = std::min(min_col_vector[pi],pj_minus);
	    max_col_vector[pi] = std::max(max_col_vector[pi],pj_plus);
	}
	
	for (size_t pi=pi_minus; pi<pi_min; pi++) {
	    max_col_vector[pi] = std::max(max_col_vector[pi],pj_max);
	}
	
	for (size_t pi=pi_max+1; pi<=pi_plus; pi++) {
	    min_col_vector[pi] = std::min(min_col_vector[pi],pj_min);
	}
	    
    }
    
#endif
    
    // print_debug(std::cout);
}

void
TraceController::TraceRange::print_debug(std::ostream & out) const {
    out << "min_col_vector: ";
    for (std::vector<size_type>::const_iterator it=min_col_vector.begin(); it!=min_col_vector.end(); ++it) {
	out.width(3);
	out << *it;
    }
    out << std::endl;
    out << "max_col_vector: ";
    for (std::vector<size_type>::const_iterator it=max_col_vector.begin(); it!=max_col_vector.end(); ++it) {
	out.width(3);
	out << *it;
    }
    out << std::endl;
}

TraceController::~TraceController() {}

void
TraceController::constrain_wo_ref(size_type lenA, size_type lenB, size_type delta) {
    // fill vectors for min_j and max_j
    for (size_type i=0; i<=lenA; i++) {
	min_col_vector[i] = std::max((size_type)0, (size_type)(ceil(i*lenB/lenA-delta)));
	max_col_vector[i] = std::min(lenB, (size_type)(floor(i*lenB/lenA+delta)));
    }
}

/* Construct from MultipleAlignment (as needed for progressive alignment) */

TraceController::TraceController(Sequence seqA, Sequence seqB, const MultipleAlignment *ma, int delta_)
  : delta(delta_) {
    
    size_type lenA = seqA.length();
    size_type lenB = seqB.length();
    
    min_col_vector.resize(lenA+1);
    max_col_vector.resize(lenA+1);
    
    // initialize vectors least constrained
    fill(min_col_vector.begin(), min_col_vector.end(),0);
    fill(max_col_vector.begin(), max_col_vector.end(),lenB);
    
    if ( delta_ == -1 ) { // no constraints!	
	// do nothing
    }
    else if ( ma == NULL ) { 
	// constraints due to delta but no reference alignment
	constrain_wo_ref(lenA, lenB, (size_type)delta);
    } else {
	
	// HERE: delta >= 0 and reference alignment ma is given
	
	
	// ----------------------------------------
	// Compute valid trace cells from REFERENCE ALIGNMENT
	//
	// constrain the valid traces from the traces
	// of all pairwise alignments between seqA and seqB
	// as given in the reference alignment
	// 

	// construct multiple alignment objects out of sequence objects
	// since this allows easier access and provides mappings pos_to_col,
	// col_to_pos
	MultipleAlignment maSeqA(seqA);
	MultipleAlignment maSeqB(seqB);
    
	//  iterate over all pairs of rows in the multiple alignment of seqA and seqB
	for (size_type i=0; i<maSeqA.size(); ++i) {
	    const MultipleAlignment::SeqEntry &seqentryA = maSeqA.seqentry(i);
	    // get alignment string in reference corresponding to seqentryA
	    const std::string &nameA = seqentryA.name();
	    const MultipleAlignment::SeqEntry &ref_seqentryA = ma->seqentry(nameA);
	
	    for (size_type j=0; j<maSeqB.size(); ++j) {
		const MultipleAlignment::SeqEntry &seqentryB = maSeqB.seqentry(j);
		// get alignment string in reference corresponding to seqentryB
		const std::string &nameB = seqentryB.name();
		const MultipleAlignment::SeqEntry &ref_seqentryB = ma->seqentry(nameB);
	    
		// construct trace for current sequences A and B
		TraceRange tr(seqentryA,seqentryB,ref_seqentryA,ref_seqentryB,delta);
		
		//std::cout << nameA << " " << nameB << std::endl;
		//tr.print_debug(std::cout);
	
		// combine existing trace range with new trace +/- delta
		merge_in_trace_range(tr);
	    
	    }
	}

#ifndef NDEBGUG
	//std::cout << "Merged:" << std::endl;
	//TraceRange::print_debug(std::cout);
	
	for (size_type i=1; i < min_col_vector.size(); ++i) {
	    assert(min_col_vector[i-1]<=min_col_vector[i]);
	    assert(max_col_vector[i-1]<=max_col_vector[i]);
	    assert(min_col_vector[i]<=max_col_vector[i]); // otherwise trace range inconsistent
	    assert(max_col_vector[i-1]+1>=min_col_vector[i]); // ranges connected/overlap, otherwise trace is inconsistent
	}
#endif
    }
    
    //TraceRange::print_debug(std::cout);
}

void
TraceController::merge_in_trace_range(const TraceRange &tr) {
    // intersect trace range of *this with trace
    for (size_type i=0; i<=tr.rows(); i++) {
	
	min_col_vector[i] = std::max( min_col_vector[i], tr.min_col(i) );
	max_col_vector[i] = std::min( max_col_vector[i], tr.max_col(i) );
	
	
	// intersecting may lead to inconsistency, check this here.
	// probably it will be necessary to replace the intersection idea
	// by a more relaxed merging strategy
	if ( min_col_vector[i] > max_col_vector[i] 
	     || 
	     ((i>0) &&  (max_col_vector[i-1]+1<min_col_vector[i]))) {
	    std::cerr << "Inconsistent trace range due to max-diff heuristic" << std::endl;
	    exit(-1); // ATTENTION: think later what to do about that
	}
    }
}
