c     include this file when failur IS used in element/material routines,
c     or in other words when a loop is multi threaded
      integer ifail
      logical failur,failgl,failst
      common/failcm/failst,failgl
C_TASKCOMMON (failcmloc)
      common/failcmloc/ifail(nlq),failur
