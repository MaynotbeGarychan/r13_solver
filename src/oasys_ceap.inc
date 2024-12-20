!
!   If any parameters in this file are modified, then please also
!   modify them in oasys_ceap_kw.inc.
!
      character*256 fn_damp
      common /ceap_damp_1/fn_damp
      integer iceap_ctl_len
      parameter (iceap_ctl_len=140)
      real ceap_r_1
      common /ceap_damp_2/  ceap_r_1(iceap_ctl_len)
      integer  ceap_i_cnt
      common /ceap_damp_3/  ceap_i_cnt(iceap_ctl_len)
!
!   The contents of ceap_i_cnt and ceap_r_1 are described in 
!   oasys_ceap_array_contents.inc. These inc files have been
!   split to prevent the need to recompile many overlays 
!   just because a comment line has changed in oasys_ceap.inc.
!
      integer iceap_on, icpf, icpf_atype, icpf_output, ipwptime, icpa,  &
     & icpc,icpa_prt
      common /iceap_ctrl/                                               &
     &                       iceap_on,                                  &
     &                       icpf,                                      &
     &                       icpf_atype,                                &
     &                       icpf_output,                               &
     &                       ipwptime,                                  &
     &                       icpa,                                      &
     &                       icpc,                                      &
     &                       icpa_prt
!
      real cpf_datum, cpf_wtable, cpf_pfrho, cpf_grav, cpf_pfbkm,       &
     &     cpf_tmf, event_time, event_mult, cpf_parho, cpf_pair,        &
     &     cpf_rta
      common /rceap_ctrl/                                               &
     &                       cpf_datum,                                 &
     &                       cpf_wtable,                                &
     &                       cpf_pfrho,                                 &
     &                       cpf_grav,                                  &
     &                       cpf_pfbkm,                                 &
     &                       cpf_tmf,                                   &
     &                       event_time,                                &
     &                       event_mult,                                &
     &                       cpf_parho,cpf_pair,cpf_rta
!
! icpf       - pore water pressure switch - 0 => off, 1 => on
! icpf_atype - pwp type - 0 => none, 
!                       - 1 => undrained,
!                       - 2 => drained,
!                       - 3 => TDC,
!                       - 4 => Consolidate to steady state,
!                       - 5 => Drained in DR, undrained in transient
! cpf_datum  - z-datum for pwp calculations
! cpf_wtable - default z coord of water table
! cpf_pfrho  - default pore fluid density
! cpf_grav   - gravity loading applied to pore fluid
! cpf_pfbkm  - default pore fluid bulk modulus
! cpf_parho  - default pore air density
! cpf_rair   - default pore air gas constant (Joule/Kg/Deg)
! cpf_tair   - default pore air temperature
! cpf_pair   - atmospheric pressure
! icpf_output- switch output to effective stress tensor
!
!
! number of pore water pressure history variables per solid element
! (now reduced from 15 to 5) - npwpcon initialized to 5 in init_once
      integer npwpcon
!      parameter (npwpcon=5)
      common/ceap_pwp1/npwpcon,npwpnod
! number of pore water pressure history variables per node          
      integer npwpnod
!      parameter (npwpnod=25)
!
!
! CEAP Arrays
!
      common /ceap_pnt/mp_cer,mp_cei
      common /ceap_ar1/ piceap,prceap
      pointer(mp_cei,iceap)
      pointer(mp_cer,rceap)
      integer iicpnt_len
      parameter (iicpnt_len=100)
      integer ircpnt_len
      parameter (ircpnt_len=100)
      integer piceap(iicpnt_len),prceap(ircpnt_len)
      integer iceap(1)
      real    rceap(1)
      integer ihcpnt_len
      parameter (ihcpnt_len=500)
      integer ihceap(ihcpnt_len)
      common/ceap_handles/ihceap
!
!     iceap      - main integer ceap array
!     rceap      - main real ceap array
!
!     mp_cei     - memory pointer to base of iceap
!     mp_cer     - memory pointer to base of rceap
!
!     iicpnt_len - length of ceap 'pointers to integers' array 
!                  (piceap)
!     ircpnt_len - length of ceap 'pointers to reals' array (prceap)
!     ihcpnt_len - length of ceap handles (ihceap)
!
!     piceap     - pointers to iceap (listed in oasys_ceap_array_contents.inc)
!     prceap     - pointers to rceap (listed in oasys_ceap_array_contents.inc)
!     ihceap     - handles to ceap storage in rmem                            
