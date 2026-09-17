program test_legacy_riccati_trajectory_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t
  use legacy_riccati_trajectory_2d, only : propagate_legacy_riccati_to_target
  implicit none

  call test_frozen_radial_source_checkpoints()

  print '(a)', "legacy Riccati 2D trajectory tests passed"

contains

  subroutine test_frozen_radial_source_checkpoints()
    type(cartesian_mesh_2d_t) :: mesh
    type(straight_trajectory_2d_t) :: trajectory
    type(trajectory_self_energy_2d_t) :: self_energy
    complex(rk), allocatable :: from_entry(:, :), from_exit(:, :)
    complex(rk) :: expected_forward(0:3, 3), expected_reverse(0:3, 3)
    integer :: point, source_index

    call make_uniform_cartesian_mesh(-0.4_rk, 0.4_rk, 4, &
                                     -0.1_rk, 0.1_rk, 1, mesh)
    call build_straight_trajectory_2d(mesh, [0.0_rk, 0.0_rk], &
                                      [1.0_rk, 0.0_rk, 0.0_rk], &
                                      0.2_rk, trajectory)
    call require(trajectory%sample_count() == 5 .and. &
                 trajectory%target_sample == 3, &
                 "frozen Riccati checkpoint trajectory has the wrong grid")

    allocate(self_energy%triplet_pair_potential(3, 5))
    allocate(self_energy%diagonal_shift(5))
    do point = 1, 5
      source_index = point - 3
      self_energy%triplet_pair_potential(1, point) = cmplx( &
        0.23_rk + 0.015_rk * real(source_index, rk), &
        0.04_rk - 0.012_rk * real(source_index, rk), kind=rk)
      self_energy%triplet_pair_potential(2, point) = cmplx( &
        -0.08_rk + 0.009_rk * real(source_index, rk), &
        0.17_rk + 0.006_rk * real(source_index, rk), kind=rk)
      self_energy%triplet_pair_potential(3, point) = cmplx( &
        0.11_rk - 0.007_rk * real(source_index, rk), &
        -0.13_rk + 0.01_rk * real(source_index, rk), kind=rk)
      self_energy%diagonal_shift(point) = &
        0.018_rk * real(source_index, rk)
    end do

    expected_forward = cmplx(0.0_rk, 0.0_rk, kind=rk)
    expected_forward(1:3, 1) = [ &
      cmplx(-4.16287674828011e-2_rk, 1.14113655080881e-1_rk, kind=rk), &
      cmplx(-8.90497686716779e-2_rk, -6.03064913261076e-2_rk, kind=rk), &
      cmplx(8.38229851553367e-2_rk, 7.51486865125393e-2_rk, kind=rk)]
    expected_forward(1:3, 2) = [ &
      cmplx(-3.90906259133229e-2_rk, 1.14416020940484e-1_rk, kind=rk), &
      cmplx(-8.89008683058641e-2_rk, -5.89276884832541e-2_rk, kind=rk), &
      cmplx(8.24471215602192e-2_rk, 7.37885530982586e-2_rk, kind=rk)]
    expected_forward(1:3, 3) = [ &
      cmplx(-3.45416635611957e-2_rk, 1.17489705712979e-1_rk, kind=rk), &
      cmplx(-9.01737612952456e-2_rk, -5.58797007439699e-2_rk, kind=rk), &
      cmplx(8.02906423983453e-2_rk, 7.10955985562741e-2_rk, kind=rk)]

    ! Reverse checkpoints are ordered by propagation: source indices 2, 1, 0.
    expected_reverse = cmplx(0.0_rk, 0.0_rk, kind=rk)
    expected_reverse(1:3, 1) = [ &
      cmplx(1.51944524958707e-2_rk, 1.49450181030876e-1_rk, kind=rk), &
      cmplx(1.03447724849656e-1_rk, -3.99082230021858e-2_rk, kind=rk), &
      cmplx(-6.11783682746152e-2_rk, 5.78452048504464e-2_rk, kind=rk)]
    expected_reverse(1:3, 2) = [ &
      cmplx(1.50939112061459e-2_rk, 1.46886364861708e-1_rk, kind=rk), &
      cmplx(1.01564261799772e-1_rk, -4.03664855254483e-2_rk, kind=rk), &
      cmplx(-6.14694115981775e-2_rk, 5.80083799869746e-2_rk, kind=rk)]
    expected_reverse(1:3, 3) = [ &
      cmplx(1.62659859791399e-2_rk, 1.42678313138293e-1_rk, kind=rk), &
      cmplx(9.95709656830779e-2_rk, -4.16922195641589e-2_rk, kind=rk), &
      cmplx(-6.37043298757304e-2_rk, 5.90495834745839e-2_rk, kind=rk)]

    call propagate_legacy_riccati_to_target( &
      trajectory, self_energy, cmplx(0.0_rk, 0.83_rk, kind=rk), 3, 0.0_rk, &
      from_entry, from_exit)

    call require(maxval(abs(from_entry(:, 1:3) - expected_forward)) < &
                 8.0e-15_rk, &
                 "forward 2D Riccati checkpoints differ from new_src")
    call require(maxval(abs(from_exit(:, [5, 4, 3]) - expected_reverse)) < &
                 8.0e-15_rk, &
                 "reverse 2D Riccati checkpoints differ from new_src")
  end subroutine test_frozen_radial_source_checkpoints


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_legacy_riccati_trajectory_2d
