SUBROUTINE calc_q3_cluster(n_atoms,n_q3_neigh,q3_cutoff,n_bonds,box,x,largest_cluster)

  IMPLICIT NONE

  ! Passed
  INTEGER,INTENT(IN) :: n_atoms
  INTEGER,INTENT(IN) :: n_q3_neigh
  INTEGER,INTENT(IN) :: n_bonds
  DOUBLE PRECISION,INTENT(IN) :: q3_cutoff
  DOUBLE PRECISION,INTENT(IN) :: box
  DOUBLE PRECISION,INTENT(IN),DIMENSION(0:3*n_atoms-1) :: x
  INTEGER,INTENT(OUT) :: largest_cluster

  ! Local
  INTEGER :: i,j,k
  INTEGER :: iatom,jatom, iindx,jindx
  DOUBLE PRECISION,DIMENSION(3) :: xi,xj,dxij
  DOUBLE PRECISION :: rij,rijsq

  ! Neighbor list
  INTEGER,DIMENSION(n_atoms) :: n_nneigh
  INTEGER,DIMENSION(75,n_atoms) :: nneigh_lst
  DOUBLE PRECISION,DIMENSION(4,75,n_atoms) :: rij_lst
  INTEGER :: nneigh_temp
  DOUBLE PRECISION,DIMENSION(4) :: rij_temp
  DOUBLE PRECISION :: q3_cutoffsq
 
  ! Spherical Harmonics
  DOUBLE PRECISION :: Pi,factor, fcij
  DOUBLE PRECISION :: pre0,pre1,pre2,pre3
  DOUBLE PRECISION :: pre(-3:3)
  DOUBLE PRECISION :: cost(3),sint(3)
  DOUBLE PRECISION :: cosp(3), sinp(3)
  COMPLEX*16 :: Y(-3:3),comaux
  COMPLEX*16 :: sum_Y(-3:3),sum_q3(-3:3,n_atoms)
  DOUBLE PRECISION,DIMENSION(75,n_atoms) :: dij
  DOUBLE PRECISION,DIMENSION(n_atoms) :: q3iatom

  ! Clustering
  INTEGER,DIMENSION(n_atoms) :: n_dij_bonds
  INTEGER :: n_clusters,n_ice
  INTEGER :: clusterID,clusterID_old, d6ijcounter
  INTEGER,DIMENSION(n_atoms) :: cluster_size, cluster_ID
  INTEGER,DIMENSION(n_atoms, n_atoms) :: cluster_lst
  LOGICAL,DIMENSION(n_atoms) :: InCluster, visited

  q3_cutoffsq =  q3_cutoff*q3_cutoff

  Pi = 3.1415926535897930d0
  factor = 2.0d0*dSQRT(Pi/7.0d0)
  fcij = factor*factor
  pre(-3) = (1.0d0/8.0d0)*dSQRT(35.0d0/Pi)
  pre(-2) = (1.0d0/4.0d0)*dSQRT(105.0d0/(2.*Pi))
  pre(-1) = (1.0d0/8.0d0)*dSQRT(21.0d0/Pi)
  pre(0) =  (1.0d0/4.0d0)*dSQRT(7.0d0/Pi)
  pre(1) = -pre(-1)
  pre(2) =  pre(-2)
  pre(3) = -pre(-3)

  ! Initialize arrays
  n_nneigh = 0
  nneigh_lst = 0
  rij_lst = 0.0d0
  dij = 0.0d0
  q3iatom = 0.0d0
  sum_q3 = 0.0d0
  dij = 0.0d0

  ! Loop over iatom
  DO iatom = 1,n_atoms-1
    iindx = 3*(iatom-1)

    ! Store iatom's coordinates
    xi(1) = x(iindx)
    xi(2) = x(iindx+1)
    xi(3) = x(iindx+2)

    ! Loop over jatom
    DO jatom = iatom+1,n_atoms
      jindx = 3*(jatom-1)

      ! Store jatom's coordinates
      xj(1) = x(jindx)
      xj(2) = x(jindx+1)
      xj(3) = x(jindx+2)

      ! Compute ij vector & separation distance
      dxij(:) = xj(:) - xi(:)
      dxij(:) = dxij(:) - dNINT(dxij(:)/box)*box
      rijsq = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)

      ! If they are within the cutoff distance
      IF(rijsq .GT. q3_cutoffsq) CYCLE

        rij = dSQRT(rijsq)

        ! Update information for iatom

        ! Update the number of neighbors for iatom
        n_nneigh(iatom) = n_nneigh(iatom) + 1

        ! Include jatom as a neighbor of iatom
        nneigh_lst(n_nneigh(iatom),iatom)  = jatom

        ! Store the distance components for iatom
        rij_lst(1, n_nneigh(iatom), iatom) = dxij(1)
        rij_lst(2, n_nneigh(iatom), iatom) = dxij(2)
        rij_lst(3, n_nneigh(iatom), iatom) = dxij(3)
        rij_lst(4, n_nneigh(iatom), iatom) = rij

        ! Update information for jatom

        ! Update the number of neighbors for jatom
        n_nneigh(jatom) = n_nneigh(jatom) + 1

        ! Include iatom as a neighbor of jatom
        nneigh_lst(n_nneigh(jatom),jatom) = iatom

        ! Store the distance components for jatom
        rij_lst(1, n_nneigh(jatom), jatom) = -dxij(1)
        rij_lst(2, n_nneigh(jatom), jatom) = -dxij(2)
        rij_lst(3, n_nneigh(jatom), jatom) = -dxij(3)
        rij_lst(4, n_nneigh(jatom), jatom) = rij


      !ENDIF

    ENDDO
  ENDDO 

  ! Check the minimum number of neighbors are present
  IF(MINVAL(n_nneigh(:)) .LT. n_q3_neigh .AND. n_q3_neigh .GT. 0) THEN
    WRITE(*,'(A,I5,A,I3,A)') 'FATAL ERROR: ', MINLOC(n_nneigh(:)), ' HAS LESS THAN ', n_q3_neigh,' NEAREST NEIGHBORS IN SUBROUTINE Q3, ADJUST Q3_CUTOFF'
    STOP
  ENDIF

  ! Check that the maximum number of neighbors is not exceeded
  IF(MAXVAL(n_nneigh(:)) .GT. 75) THEN
    WRITE(*,'(A,I5,A)') 'FATAL ERROR: ', MAXLOC(n_nneigh(:)), ' HAS MORE THAN 75 NEAREST NEIGHBORS IN SUBROUTINE Q3'
    STOP
  ENDIF

  ! Sort the neighbors using insertion sort if using only the closest n_q3_neigh
  ! neighbors
  IF(n_q3_neigh .GT. 0 ) THEN
    DO iatom = 1,n_atoms
      DO k = 2,n_nneigh(iatom)
        nneigh_temp = nneigh_lst(k,iatom)
        rij_temp(:) = rij_lst(:,k,iatom)
        DO j = k-1,1,-1
          IF(rij_lst(4,j,iatom) .LE. rij_temp(4)) GOTO 10
          nneigh_lst(j + 1, iatom) = nneigh_lst(j,iatom)
          rij_lst(:,j + 1,iatom) = rij_lst(:, j, iatom)
        ENDDO
        jatom = 0
        10 CONTINUE
        nneigh_lst(j + 1, iatom) = nneigh_temp
        rij_lst(:,j + 1, iatom) = rij_temp(:)
      ENDDO
      n_nneigh(iatom) = n_q3_neigh
    ENDDO
  ENDIF


  ! Calculate q3 and dij

  ! Loop over the atoms
  DO iatom = 1,n_atoms
    iindx = 3*(iatom-1)

    ! Zero harmonic sum
    sum_Y = 0.0d0

    ! Loop over the nearest neighbors of iatom
    DO j = 1, n_nneigh(iatom)
      jatom = nneigh_lst(j,iatom)
      jindx = 3*(jatom-1)

      ! Get separation vector for list
      dxij(1) = rij_lst(1,j,iatom)
      dxij(2) = rij_lst(2,j,iatom)
      dxij(3) = rij_lst(3,j,iatom)
      rij = rij_lst(4,j,iatom)

      ! Calculate cos(theta) and sin(theta)
      cost(1) = dxij(3)/rij
      sint(1) = (1.0d0-cost(1)*cost(1))**0.5d0

      ! Calculate cos(phi) and sin(phi)
      IF (dABS(sint(1)) .LT. 1.0d-9) THEN
        cosp(1) = 1.0d0
        sinp(1) = 0.0d0
      ELSE
        cosp(1) = dxij(1)/(rij*sint(1))
        sinp(1) = dxij(2)/(rij*sint(1))
      ENDIF

      ! Calculate powers of sin(theta)
      sint(2) = sint(1)*sint(1)
      sint(3) = sint(2)*sint(1)

      ! Calculate powers of cos(theta)
      cost(2) = cost(1)*cost(1)
      cost(3) = cost(2)*cost(1)

      ! Calculate sin/cos (2phi)
      sinp(2) = 2.0d0*sinp(1)*cosp(1)
      cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

      ! Calculate sin/cos (3phi)
      sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
      cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)

      ! Update the harmonic sum
      Y(-3) = pre(-3)*sint(3)*dCMPLX(cosp(3),-sinp(3))
      Y(-2) = pre(-2)*sint(2)*cost(1)*dCMPLX(cosp(2),-sinp(2))
      Y(-1) = pre(-1)*sint(1)*(5.0d0*cost(2)-1.0)*dCMPLX(cosp(1),-sinp(1))
      Y(0) = pre(0)*(5.*cost(3)-3.0d0*cost(1))
      Y(1) = pre(1)*sint(1)*(5.0d0*cost(2)-1.0)*dCMPLX(cosp(1),sinp(1))
      Y(2) = pre(2)*sint(2)*cost(1)*dCMPLX(cosp(2),sinp(2))
      Y(3) = pre(3)*sint(3)*dCMPLX(cosp(3),sinp(3))

      sum_Y = sum_Y + Y

    ENDDO

    sum_q3(:,iatom) = sum_Y(:)

    ! Calculate q3 for each atom
    DO i = -3,3
      q3iatom(iatom) = q3iatom(iatom) + sum_q3(i,iatom)*dCONJG(sum_q3(i,iatom))
    ENDDO
    q3iatom(iatom) = factor*dSQRT(q3iatom(iatom))/DBLE(n_nneigh(iatom))
  ENDDO

  ! Loop over the atoms
  n_dij_bonds = 0
  n_ice = 0
  DO iatom = 1,n_atoms
    DO j = 1, n_nneigh(iatom)
      jatom = nneigh_lst(j,iatom)
      DO i = -3,3
        dij(j,iatom) = dij(j,iatom) + sum_q3(i,iatom)*CONJG(sum_q3(i,jatom))
      ENDDO
      dij(j,iatom) = fcij*dij(j,iatom)/(q3iatom(iatom)*q3iatom(jatom))/DBLE(n_nneigh(iatom)*n_nneigh(jatom))
      IF((dij(j,iatom) .LE. -0.82d0) .OR. ((dij(j,iatom) .GE. -0.145d0) .AND. (dij(j,iatom) .LE. -0.065d0))) THEN
         n_dij_bonds(iatom) = n_dij_bonds(iatom) + 1
      ENDIF
    ENDDO
  ENDDO


  ! Initialize variables
  InCluster = .FALSE. !FALSE if atom i not in a cluster 
  n_clusters = 0
  cluster_size = 0
  cluster_lst = 0
  cluster_ID = 0

  !Loop over the atoms
  DO iatom = 1,n_atoms

    clusterID = 0

    IF (n_dij_bonds(iatom) .LT. n_bonds) CYCLE

      InCluster(iatom) = .TRUE.
  
      !Search iatoms's n.n. to see if they are in cluster(s)
      DO j = 1, n_nneigh(iatom)
        jatom = nneigh_lst(j,iatom)
 
        IF (InCluster(jatom) .EQV. .TRUE.) THEN !one of iatom's n.n. is in a cluster

          IF (clusterID .NE. 0) THEN !another one of iatoms's n.n. is also in a cluster

            IF (clusterID .NE. cluster_ID(jatom)) THEN !two of iatoms's n.n. belong to "different" clusters...combine

              !keep lower cluster ID number
              IF (clusterID .LT. cluster_ID(jatom)) THEN
                clusterID_old = cluster_ID(jatom)
              ELSE
                clusterID_old = clusterID
                clusterID = cluster_ID(jatom)
              ENDIF

              ! Merge clusters, keeping lower cluster ID number
              DO k = 1,cluster_size(clusterID_old)
                cluster_lst(cluster_size(clusterID) + k,clusterID) = cluster_lst(k,clusterID_old)
                cluster_ID(cluster_lst(k,clusterID_old)) = clusterID
              ENDDO
              cluster_size(clusterID) = cluster_size(clusterID) + cluster_size(clusterID_old)
              cluster_lst(:,clusterID_old) = 0 ! Clear old cluster data, reduce count
              cluster_size(clusterID_old) = 0     
              n_clusters = n_clusters - 1      

              !Place LAST cluster in cluster list into newly empty slot
              IF (clusterID_old .NE. (n_clusters+1)) THEN ! Must do move
                DO k = 1,cluster_size(n_clusters+1)
                 cluster_lst(k,clusterID_old) = cluster_lst(k,(n_clusters+1))
                 cluster_ID(cluster_lst(k,clusterID_old)) = clusterID_old
                ENDDO
                cluster_size(clusterID_old) = cluster_size(n_clusters+1)  ! Clear data for moved cluster 
                cluster_lst(:,(n_clusters+1)) = 0 
                cluster_size(n_clusters+1) = 0    
              ENDIF
            ENDIF
       
          ELSE
            clusterID = cluster_ID(jatom)
          ENDIF

        ENDIF
      ENDDO  ! j neighbor

      IF (clusterID .EQ. 0) then !if none of iatom's n.n. are in a cluster, imol begins new cluster
        n_clusters = n_clusters + 1
        clusterID = n_clusters
      ENDIF

      !Add iatom to the cluster
      cluster_ID(iatom) = clusterID
      cluster_size(clusterID) = cluster_size(clusterID) + 1
      cluster_lst(cluster_size(clusterID),clusterID) = iatom !cluster_lst(cluster number, list of atoms in cluster)

  ENDDO ! iatom
  
  largest_cluster = MAXVAL(cluster_size)
  
END SUBROUTINE
