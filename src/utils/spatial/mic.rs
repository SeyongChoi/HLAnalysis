pub fn mic(
    vector:[f64;3],
    cell: &[[f64;3];3],
    cell_inv: &[[f64;3];3]
)->[f64;3]{
    // Fractional 좌표로 변환
    let fractional = [
        cell_inv[0][0] * vector[0] + cell_inv[1][0] * vector[1] + cell_inv[2][0] * vector[2],
        cell_inv[0][1] * vector[0] + cell_inv[1][1] * vector[1] + cell_inv[2][1] * vector[2],
        cell_inv[0][2] * vector[0] + cell_inv[1][2] * vector[1] + cell_inv[2][2] * vector[2],
    ];

    // Fractional 좌표를 -0.5 ~ 0.5 범위로 wrap
    let fractional_wrapped = fractional.map(|x| x - (x + 0.5).floor());

    // 다시 real space로 변환
    [
        cell[0][0] * fractional_wrapped[0] + cell[1][0] * fractional_wrapped[1] + cell[2][0] * fractional_wrapped[2],
        cell[0][1] * fractional_wrapped[0] + cell[1][1] * fractional_wrapped[1] + cell[2][1] * fractional_wrapped[2],
        cell[0][2] * fractional_wrapped[0] + cell[1][2] * fractional_wrapped[1] + cell[2][2] * fractional_wrapped[2],
    ]
}
