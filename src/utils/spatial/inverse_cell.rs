pub fn inverse_cell(cell: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    // Step 1: 행렬식 계산
    let det = cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
        - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
        + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0]);

    if det == 0.0 {
        panic!("Matrix is singular and cannot be inverted"); // 역행렬이 존재하지 않으면 프로그램 종료
    }

    // Step 2: 여인수 행렬 계산
    let cofactor = [
        [
            cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1],
            -(cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0]),
            cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0],
        ],
        [
            -(cell[0][1] * cell[2][2] - cell[0][2] * cell[2][1]),
            cell[0][0] * cell[2][2] - cell[0][2] * cell[2][0],
            -(cell[0][0] * cell[2][1] - cell[0][1] * cell[2][0]),
        ],
        [
            cell[0][1] * cell[1][2] - cell[0][2] * cell[1][1],
            -(cell[0][0] * cell[1][2] - cell[0][2] * cell[1][0]),
            cell[0][0] * cell[1][1] - cell[0][1] * cell[1][0],
        ],
    ];

    // Step 3: 전치 행렬 (Adjugate) 계산
    let adjugate = [
        [cofactor[0][0], cofactor[1][0], cofactor[2][0]],
        [cofactor[0][1], cofactor[1][1], cofactor[2][1]],
        [cofactor[0][2], cofactor[1][2], cofactor[2][2]],
    ];

    // Step 4: 역행렬 계산
    let inv_det = 1.0 / det;
    let mut inv_cell = [[0.0; 3]; 3];

    for i in 0..3 {
        for j in 0..3 {
            inv_cell[i][j] = adjugate[i][j] * inv_det;
        }
    }

    inv_cell
}