// Copyright (c) 2026 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_MARK_MANAGEMENT_H
#define CGAL_MARK_MANAGEMENT_H 1

#include <CGAL/assertions.h>
#include <bitset>

namespace CGAL {
  /////////////////////////////////////////////////////////////////////////////
  template<typename Refs, typename Refs::size_type NUMBER_OF_MARKS=32>
  class Mark_management_basic
  {
  public:
    class Exception_no_more_available_mark {};

    using size_type=typename Refs::size_type;

    /// Number of marks
    static const size_type NB_MARKS=NUMBER_OF_MARKS;

    Mark_management_basic()
    { reset(); }

    // Reset all the data used to manage marks
    void reset()
    {
      mnb_used_marks=0;
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        mfree_marks_stack[i]   =i;
        mindex_marks[i]        =i;
        mnb_marked_darts[i]    =0;
        mnb_times_used_marks[i]=0;
      }
    }

    // Call after a clear of all the darts: only reset the number of marked
    // darts, but keep the used marks as it.
    void clear_darts()
    {
      for ( size_type i=0; i<NB_MARKS; ++i)
      { mnb_marked_darts[i]=0; }
    }

    template<typename Mark_management2>
    void copy(const Mark_management2& other)
    {
      // Duplicate all marks of amap not yet used
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        if(!is_used(i) && other.is_used(i))
        {
          CGAL_assertion(mnb_used_marks<NB_MARKS);
          // 1) Remove mark i from mfree_marks_stack
          //     (replace it by the last free mark)
          mfree_marks_stack[mindex_marks[i]]=
              mfree_marks_stack[NB_MARKS-mnb_used_marks-1];
          mindex_marks[mfree_marks_stack[mindex_marks[i]]]=mindex_marks[i];
          // 2) Update use mark stack
          mused_marks_stack[mnb_used_marks]=i;
          mindex_marks[i]=mnb_used_marks;
          mnb_times_used_marks[i]=1;
          ++mnb_used_marks;
        }
      }
    }

    void swap(Mark_management_basic& other)
    {
      std::swap(mnb_times_used_marks, other.mnb_times_used_marks);
      std::swap(mnb_used_marks, other.mnb_used_marks);
      std::swap(mindex_marks, other.mindex_marks);
      std::swap(mfree_marks_stack, other.mfree_marks_stack);
      std::swap(mused_marks_stack, other.mused_marks_stack);
      std::swap(mnb_marked_darts, other.mnb_marked_darts);
    }

    /** Count the number of used marks.
     * @return the number of used marks.
     */
    size_type number_of_used_marks() const
    { return mnb_used_marks; }

    /** Tests if a given mark is used.
     *  @return true iff the mark is used.
     */
    bool is_used(size_type amark) const
    {
      CGAL_assertion(amark<NB_MARKS);
      return (mnb_times_used_marks[amark]!=0);
    }

    /** Get a new free mark and return its index.
     * All the darts are unmarked for this mark.
     * @return the index of the new mark.
     * @pre mnb_used_marks < NB_MARKS
     */
    size_type get_new_mark(std::size_t /*nb_darts*/) const
    {
      if (mnb_used_marks==NB_MARKS)
      {
        std::cerr << "Not enough Boolean marks: "
                     "increase NB_MARKS in item class." << std::endl;
        std::cerr << "  (exception launched)" << std::endl;
        throw Exception_no_more_available_mark();
      }

      size_type m=mfree_marks_stack[mnb_used_marks];
      mused_marks_stack[mnb_used_marks]=m;

      mindex_marks[m]=mnb_used_marks;
      mnb_times_used_marks[m]=1;

      ++mnb_used_marks;

      return m;
    }

    /** Increase the number of times a mark is used.
     *  @param amark the mark to share.
     */
    void share_a_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      ++mnb_times_used_marks[amark];
    }

    /** @return the number of times a mark is used.
     *  @param amark the mark to share.
     */
    size_type get_number_of_times_mark_used(size_type amark) const
    {
      CGAL_assertion( amark<NB_MARKS );
      return mnb_times_used_marks[amark];
    }

    /** Negate the mark of all the darts for a given mark.
     * After this call, all the marked darts become unmarked and all the
     * unmarked darts become marked (in constant time operation).
     * @param amark the mark index
     */
    void negate_mark(size_type amark, std::size_t nb_darts) const
    {
      CGAL_assertion( is_used(amark) );
      mnb_marked_darts[amark]=nb_darts-mnb_marked_darts[amark];
    }

    /** Free a given mark.
     * @param amark the given mark.
     */
    void free_mark(size_type amark) const
    {
      CGAL_assertion( is_used(amark) );

      if ( mnb_times_used_marks[amark]>1 )
      {
        --mnb_times_used_marks[amark];
        return;
      }

      // 1) We remove amark from the array mused_marks_stack by
      //    replacing it with the last mark in this array.
      mused_marks_stack[mindex_marks[amark]] =
          mused_marks_stack[--mnb_used_marks];
      mindex_marks[mused_marks_stack[mnb_used_marks]] =
          mindex_marks[amark];

      // 2) We add amark in the array mfree_marks_stack and update its index.
      mfree_marks_stack[ mnb_used_marks ] = amark;
      mindex_marks[amark] = mnb_used_marks;

      mnb_times_used_marks[amark]=0;
    }

    /**  Count the number of marked darts for a given mark.
     * @param amark the mark index.
     * @return the number of marked darts for amark.
     */
    size_type number_of_marked_darts(size_type amark) const
    {
      CGAL_assertion(is_used(amark));
      return mnb_marked_darts[amark];
    }

      protected:
    /// Number of times each mark is used. 0 if the mark is free.
    mutable std::array<size_type, NB_MARKS> mnb_times_used_marks{0};

    /// Number of used marks.
    mutable size_type mnb_used_marks=0;

    /// Index of each mark, in mfree_marks_stack or in mfree_marks_stack.
    mutable std::array<size_type, NB_MARKS> mindex_marks{0};

    /// "Stack" of free marks.
    mutable std::array<size_type, NB_MARKS> mfree_marks_stack{0};

    /// "Stack" of used marks.
    mutable std::array<size_type, NB_MARKS> mused_marks_stack{0};

    /// Number of marked darts for each used marks.
    mutable std::array<size_type, NB_MARKS> mnb_marked_darts{0};
  };
  /////////////////////////////////////////////////////////////////////////////
  template<typename Refs, typename Refs::size_type NUMBER_OF_MARKS=32>
  class Mark_management_with_bitset: public Mark_management_basic<Refs, 32>
  {
  public:
    using Base=Mark_management_basic<Refs, 32>;
    using size_type=typename Refs::size_type;
    using Base::NB_MARKS;

    Mark_management_with_bitset()
    { reset(); }

    void reset()
    {
      Base::reset();
      mmask_marks.reset();
    }

    template<typename Mark_management2>
    void copy(const Mark_management2& other)
    {
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        if(!this->is_used(i) && other.is_used(i))
        { mmask_marks[i]=other.mmask_marks[i]; }
      }
      Base::copy(other);
    }

    void swap(Mark_management_with_bitset<Refs>& other)
    {
      Base::swap(other);
      std::swap(mmask_marks, other.mmask_marks);
    }

    void negate_mark(size_type amark, std::size_t nb_darts) const
    {
      Base::negate_mark(amark, nb_darts);
      mmask_marks.flip(amark);
    }

    /** Get the mask associated to a given mark.
     * @param amark the mark.
     * @return the mask associated to mark amark.
     */
    bool get_mask_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      return mmask_marks[amark];
    }

    /** Unmark all the darts of the map for a given mark if this is possible
     *  in constant time
     * @param amark the given mark.
     */
    bool unmark_all_if_possible(size_type amark, std::size_t nb_darts) const
    {
      CGAL_assertion(this->is_used(amark));

      if (this->mnb_marked_darts[amark]==nb_darts)
      { negate_mark(amark, nb_darts); } // all darts are marked

      return this->mnb_marked_darts[amark]==0;
    }

  protected:
    /// Mask marks to know the value of unmark dart, for each index i.
    mutable std::bitset<NB_MARKS> mmask_marks;
  };
  /////////////////////////////////////////////////////////////////////////////
  template<typename Refs>
  class Mark_management_bitset_on_dart:
           public Mark_management_with_bitset<Refs, 32>
  {
  public:
    using Base=Mark_management_with_bitset<Refs, 32>;
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;
    using Base::NB_MARKS;

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& storage, Dart_descriptor ADart) const
    { set_dart_marks(storage, ADart, this->mmask_marks); }

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart) const
    {
      // We update the number of marked darts.
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if (is_marked(storage, ADart, this->mused_marks_stack[i]))
        { --this->mnb_marked_darts[this->mused_marks_stack[i]]; }
      }
    }

    /// copy the marks of ADart1 on ADart2
    template<typename Refs2, typename Mark_management2>
    void copy_marks_of_dart(const Refs& storage, Dart_const_descriptor ADart1,
                            const Refs2& storage2, const Mark_management2& mm2,
                            typename Refs2::Dart_const_descriptor ADart2) const
    {
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if(is_marked(storage, ADart1, this->mused_marks_stack[i]))
        { mm2.set_mark_to(storage2, ADart2, this->mused_marks_stack[i], true); }
      }
    }

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& /* storage */,
                       Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& /* storage */,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    {
      CGAL_assertion(ADart!=nullptr);
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& /* storage */,
                        Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion(ADart!=nullptr);
      ADart->flip_mark(amark);
    }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const Refs& /* storage */,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      CGAL_assertion(ADart!=nullptr);
      ADart->set_marks(amarks);
    }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const Refs& /* storage */,
                                         Dart_const_descriptor ADart) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_marks();
    }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const Refs& storage,
                   Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ this->mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const Refs& storage,
                                    Dart_const_descriptor adart) const
    { return get_dart_marks(storage, adart) ^ this->mmask_marks; }

    /** Tests if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(const Refs& storage,
                   Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion(this->is_used(amark));
      return get_dart_mark(storage, adart, amark)!=this->get_mask_mark(amark);
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(const Refs& storage,
                     Dart_const_descriptor adart, size_type amark,
                     bool astate) const
    {
      CGAL_assertion(this->is_used(amark) );

      if (is_marked(storage, adart, amark)!=astate)
      {
        if (astate) { ++this->mnb_marked_darts[amark]; }
        else { --this->mnb_marked_darts[amark]; }

        flip_dart_mark(storage, adart, amark);
      }
    }

  };
  /////////////////////////////////////////////////////////////////////////////
  template<typename Refs>
  class Mark_management_bitset_on_dart_with_index:
             public Mark_management_with_bitset<Refs, 32>
  {
  public:
    using Base=Mark_management_with_bitset<Refs, 32>;
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;
    using Base::NB_MARKS;

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& storage, Dart_descriptor ADart) const
    { set_dart_marks(storage, ADart, this->mmask_marks); }

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart) const
    {
      // We update the number of marked darts.
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if (is_marked(storage, ADart, this->mused_marks_stack[i]))
        { --this->mnb_marked_darts[this->mused_marks_stack[i]]; }
      }
    }

    /// copy the marks of ADart1 on ADart2
    template<typename Refs2, typename Mark_management2>
    void copy_marks_of_dart(const Refs& storage, Dart_const_descriptor ADart1,
                            const Refs2& storage2, const Mark_management2& mm2,
                            typename Refs2::Dart_const_descriptor ADart2) const
    {
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if(is_marked(storage, ADart1, this->mused_marks_stack[i]))
        { mm2.set_mark_to(storage2, ADart2, this->mused_marks_stack[i], true); }
      }
    }

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    { return storage.mdarts[ADart].get_mark(amark); }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    { storage.mdarts[ADart].set_mark(amark, avalue); }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    { storage.mdarts[ADart].flip_mark(amark); }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(const Refs& storage,
                        Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    { storage.mdarts[ADart].set_marks(amarks); }

    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(const Refs& storage,
                                         Dart_const_descriptor ADart) const
    { return storage.mdarts[ADart].get_marks(); }

    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(const Refs& storage,
                   Dart_const_descriptor adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(storage, adart, amarks ^ this->mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(const Refs& storage,
                                    Dart_const_descriptor adart) const
    { return get_dart_marks(storage, adart) ^ this->mmask_marks; }

    /** Tests if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(const Refs& storage,
                   Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion(this->is_used(amark));
      return get_dart_mark(storage, adart, amark)!=this->get_mask_mark(amark);
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(const Refs& storage,
                     Dart_const_descriptor adart, size_type amark,
                     bool astate) const
    {
      CGAL_assertion(this->is_used(amark) );

      if (is_marked(storage, adart, amark)!=astate)
      {
        if (astate) { ++this->mnb_marked_darts[amark]; }
        else { --this->mnb_marked_darts[amark]; }

        flip_dart_mark(storage, adart, amark);
      }
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template<typename Refs>
  class Mark_management_bool_vector_in_map:
       public Mark_management_with_bitset<Refs, 32>
  {
  public:
    using Base=Mark_management_with_bitset<Refs, 32>;
    using size_type=typename Refs::size_type;
    using Dart_descriptor=typename Refs::Dart_descriptor;
    using Dart_const_descriptor=typename Refs::Dart_const_descriptor;
    using Base::NB_MARKS;

    Mark_management_bool_vector_in_map()
    { reset(); }

    void reset()
    {
      m_nb_reserved_marks=0;
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        marray_of_marks[i].clear();
        marray_of_reserved_marks[i]=0;
        mis_reserved[i]=false;
      }
      Base::reset();
    }

    template<typename Mark_management2>
    void copy(const Mark_management2& other)
    { // Do not copy reserved marks, clear them: TODO better ??
      other.m_nb_reserved_marks=0;
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        if(other.mis_reserved[i])
        {
          other.marray_of_marks[i].clear();
          other.marray_of_reserved_marks[i]=0;
          static_cast<const Mark_management2::Base&>(other).free_mark(i);
          other.mis_reserved[i]=false;
        }
      }
      Base::copy(other);
    }

    void swap(Mark_management_bool_vector_in_map& other)
    {
      Base::swap(other);
      std::swap(m_nb_reserved_marks, other.m_nb_reserved_marks);
      std::swap(mis_reserved, other.mis_reserved);
      for (size_type i=0; i<NB_MARKS; ++i)
      {
        std::swap(marray_of_marks[i], other.marray_of_marks[i]);
        std::swap(marray_of_reserved_marks[i], other.marray_of_reserved_marks[i]);
      }
    }

    size_type number_of_used_marks() const
    { return Base::number_of_used_marks()-m_nb_reserved_marks; }

    /** Reserve a new mark.
     * Get a new free mark and return its index.
     * All the darts are unmarked for this mark.
     * @return the index of the new mark.
     * @pre mnb_used_marks < NB_MARKS
     */
    size_type get_new_mark(std::size_t nb_darts) const
    {
      if(m_nb_reserved_marks>0)
      {
        mis_reserved[marray_of_reserved_marks[m_nb_reserved_marks-1]]=false;
        return marray_of_reserved_marks[--m_nb_reserved_marks];
      }

      // static std::size_t nb_full_init=0;
      // std::cout<<"get_new_mark: nb full init: "<<++nb_full_init<<std::endl;
      size_type amark=Base::get_new_mark(nb_darts);
      if(marray_of_marks[amark].size()!=nb_darts+1)
      { marray_of_marks[amark].resize(nb_darts+1, this->mmask_marks[amark]); }
      marray_of_marks[amark].assign(marray_of_marks[amark].size(),
                                    this->get_mask_mark(amark));
      return amark;
    }
/*    size_type get_new_mark(std::size_t nb_darts) const
    {
      size_type amark=Base::get_new_mark(nb_darts);
      //if(marray_of_marks[amark].size()!=nb_darts+1)
      { marray_of_marks[amark].resize(nb_darts+1, this->mmask_marks[amark]); }
      marray_of_marks[amark].assign(marray_of_marks[amark].size(),
                                    this->get_mask_mark(amark));
      return amark;
    } */

    void free_mark(size_type amark) const
    {
      if (this->mnb_times_used_marks[amark]>1)
      {
        --this->mnb_times_used_marks[amark];
        return;
      }

      // static std::size_t max_reserved=0;

      if(m_nb_reserved_marks>7)
      { Base::free_mark(amark); }
      else
      {
        mis_reserved[amark]=true;
        marray_of_reserved_marks[m_nb_reserved_marks++]=amark;
        /* if(m_nb_reserved_marks>max_reserved)
        {
          max_reserved=m_nb_reserved_marks;
          std::cout<<"Max reserved: "<<max_reserved<<std::endl;
        } */
      }
    }

    /// on_new_dart is called when a new dart is created.
    void on_new_dart(const Refs& /* storage */, Dart_descriptor ADart)
    {
      // We update the number of marked darts.
      // Here we also iterate through reserved marks (which are internally
      // considered non-free)
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if(marray_of_marks[this->mused_marks_stack[i]].size()<=ADart)
        { marray_of_marks[this->mused_marks_stack[i]].resize
              (ADart+1, this->mmask_marks[this->mused_marks_stack[i]]); }
        else
        {
          marray_of_marks[this->mused_marks_stack[i]][ADart]=
              this->mmask_marks[this->mused_marks_stack[i]];
        }
      }
    }

    /// on_delete_dart is called when a dart is deleted
    void on_delete_dart(const Refs& storage, Dart_descriptor ADart)
    {
      // We update the number of marked darts.
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if (is_marked(storage, ADart, this->mused_marks_stack[i]))
        { --this->mnb_marked_darts[this->mused_marks_stack[i]]; }
      }
    }

    /// copy the marks of ADart1 on ADart2
    template<typename Refs2, typename Mark_management2>
    void copy_marks_of_dart(const Refs& storage, Dart_const_descriptor ADart1,
                            const Refs2& storage2, const Mark_management2& mm2,
                            typename Refs2::Dart_const_descriptor ADart2) const
    {
      for (size_type i=0; i<this->mnb_used_marks; ++i)
      {
        if(is_marked(storage, ADart1, this->mused_marks_stack[i]))
        { mm2.set_mark_to(storage2, ADart2, this->mused_marks_stack[i], true); }
      }
    }

    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart, size_type amark) const
    { return marray_of_marks[amark][ADart]; }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(const Refs& storage,
                       Dart_const_descriptor ADart,
                       size_type amark, bool avalue) const
    { marray_of_marks[amark][ADart]=avalue; }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(const Refs& storage,
                        Dart_const_descriptor ADart, size_type amark) const
    { marray_of_marks[amark][ADart]=!marray_of_marks[amark][ADart]; }

    /** Tests if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(const Refs& storage,
                   Dart_const_descriptor adart, size_type amark) const
    {
      CGAL_assertion(this->is_used(amark));
      return get_dart_mark(storage, adart, amark)!=this->get_mask_mark(amark);
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(const Refs& storage,
                     Dart_const_descriptor adart, size_type amark,
                     bool astate) const
    {
      CGAL_assertion(this->is_used(amark) );

      if (is_marked(storage, adart, amark)!=astate)
      {
        if (astate) { ++this->mnb_marked_darts[amark]; }
        else { --this->mnb_marked_darts[amark]; }

        flip_dart_mark(storage, adart, amark);
      }
    }

    /** Unmark all the darts of the map for a given mark if this is possible
     *  in constant time
     * @param amark the given mark.
     */
    bool unmark_all_if_possible(size_type amark, std::size_t nb_darts) const
    {
      // this->mnb_marked_darts[amark]=0;
      if(!Base::unmark_all_if_possible(amark, nb_darts))
      {
        marray_of_marks[amark].assign(marray_of_marks[amark].size(),
                                      this->get_mask_mark(amark));
        this->mnb_marked_darts[amark]=0;
      }
      return true;
    }

  protected:
    /// array of vector of marks
    mutable std::array<std::vector<bool>, NB_MARKS> marray_of_marks;
    mutable std::size_t m_nb_reserved_marks=0;
    mutable std::array<size_type, NB_MARKS> marray_of_reserved_marks{0};
    /// Bitset to know which mark is reserved.
    mutable std::bitset<NB_MARKS> mis_reserved{false};

  };

} // namespace CGAL

#endif // CGAL_MARK_MANAGEMENT_H //
// EOF //
