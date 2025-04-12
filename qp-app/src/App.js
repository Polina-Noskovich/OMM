import React, { useState } from "react";
import * as math from "mathjs";
import { 
  Container, Typography, Button, Paper, Box, Grid, TextField, 
  Table, TableBody, TableCell, TableContainer, TableHead, 
  TableRow, Stepper, Step, StepLabel, Divider 
} from "@mui/material";

// следующая блочная матрица обратима
// H = (D* -- A(b*)')
//     (A(b*) --   0)
// где D* — это подматрица матрицы D, составленная из элементов,
// стоящих на пересечении строк и столбцов с индексами из множества J(b*);
// A(b*) — матрица, состоящая из столбцов матрицы A с индексами из множества J(b*)

const isNegative = (value) => {
  for (let i = 0; i < value.length; i++) {
    if (value[i] < 0) {
      return [value[i], i];
    }
  }
  return [null, null];
};

// Извлечение подматрицы из матрицы D
const getSubmatrixD = (D, J_) => {
  return math.matrix(
    J_.map((i) => J_.map((j) => D[i][j]))
  );
};

// Получение вектора b*
const getVectorBStar = (j0, J_, D, A) => {
  const A_j = A.map((row) => row[j0]);
  const D_j0 = D.map((row) => row[j0]);
  const D_j = J_.map((i) => D_j0[i]);
  return [...D_j, ...A_j];
};

// Обновление множества Jb и Jb*

// # Множество Jb называется опорой ограничений, а множество J(b*) — расширенной опорой ограничений.

// найдется подмножество Jb* множества индексов переменных
// для каждого индекса j ∈ Jbзвездочка выполняется ∆j (x) = 0, где 
// c'(x) = c' + x'*D; u'(x) = -cb'(x)*Ab^-1; ∆'(x) = u'(x)*A + c'(x)
const update = (Jb, JbStar, j0, teta_j0, B, A) => {
  if (j0 === teta_j0) {
    JbStar.push(teta_j0);
    return;
  }

  const JbStarNotInJb = JbStar.filter((i) => !Jb.includes(i));

  if (JbStarNotInJb.includes(teta_j0)) {
    JbStar.splice(JbStar.indexOf(teta_j0), 1);
  } else if (Jb.includes(teta_j0)) {
    const s = Jb.indexOf(teta_j0);
    let existsIndex = false;
    let jPos = -1;

    for (let i of JbStarNotInJb) {
      const value = math.dot(B[s], A.map((row) => row[i]));
      if (value !== 0) {
        existsIndex = true;
        jPos = i;
        break;
      }
    }

    if (existsIndex) {
      JbStar.splice(JbStar.indexOf(teta_j0), 1);
      Jb[s] = jPos;
      return;
    }

    Jb[s] = j0;
    JbStar[JbStar.indexOf(teta_j0)] = j0;
  }
};

// Допустимый план x задачи квадратичного программирования называется правильным опорным планом, если суще
// Формирование вектора l
const createVectorL = (D, B, Ab, J_, j0, A) => {
  const l = Array(D.length).fill(0);
  l[j0] = 1;

  const submatrixD = getSubmatrixD(D, J_);
  const matrixZ = math.zeros(J_.length, J_.length);

  const H = math.concat(
    math.concat(submatrixD, B, 1),
    math.concat(Ab, matrixZ, 1),
    0
  );

  const bStar = math.matrix(getVectorBStar(j0, J_, D, A));
  const H_inv = math.inv(H);

  const xItems = math.multiply(H_inv, bStar)._data.map((val) => -val);

  let i = 0;
  for (let item of J_) {
    l[item] = xItems[i];
    i++;
  }

  return l;
};

// Подсчет тет
const getTetas = (beta, j0, delta_x, J_, l, x) => {
  const tetas = Array(delta_x.length).fill(Infinity);

  if (beta !== 0) {
    tetas[j0] = Math.abs(delta_x[j0]) / beta;
  }

  for (let j of J_) {
    if (l[j] < 0) {
      tetas[j] = (-x[j]) / l[j];
    }
  }

  return tetas;
};

// Нахождение минимальной теты и его индекса
const minTetasAndIndex = (tetas) => {
  const tetaMin = Math.min(...tetas);
  const index = tetas.indexOf(tetaMin);
  return [tetaMin, index];
};

const QuadraticProgrammingApp = () => {
  const [activeStep, setActiveStep] = useState(0);
  const [steps, setSteps] = useState([]);
  const [result, setResult] = useState(null);
  const [params, setParams] = useState({
    A: [
      [2, 1, 0],
      [1, 0, 1],
    ],
    D: [
      [6, -2, 0],
      [-2, 2, 0],
      [0, 0, 2],
    ],
    b: [2, 1],
    c: [-1, -1, -1],
    J_on: [1, 2],
    J_ast: [1, 2],
    x: [0, 2, 1],
  });

  const handleParamChange = (param, index, subIndex, value) => {
    const newParams = { ...params };
    if (subIndex !== undefined) {
      newParams[param][index][subIndex] = Number(value);
    } else {
      newParams[param][index] = Number(value);
    }
    setParams(newParams);
  };

  const handleRun = () => {
    setActiveStep(0);
    setSteps([]);
    setResult(null);

    const { A, D, b, c, J_on, J_ast, x } = params;
    const executionSteps = [];

    const quadraticProgramming = (A, D, b, c, J_on, J_ast, x) => {
      let iteration = 0;
      while (true) {
        iteration++;
        executionSteps.push({
          iteration,
          message: `Начало итерации ${iteration}`,
          x: [...x],
          J_on: [...J_on],
          J_ast: [...J_ast],
        });

        const Ab = J_on.map((index) => A.map((row) => row[index]));
        const B = math.inv(math.matrix(Ab));
        const c_x = c.map((val, i) => val + math.dot(D[i], x));
        const c_b = J_on.map((j) => c_x[j]);
        const u = math.multiply(-1, math.multiply(c_b, B));
        const delta_x = math.add(math.multiply(u, A), c_x);

        executionSteps.push({
          iteration,
          message: "Вычислен delta_x",
          delta_x: delta_x._data ? delta_x._data : delta_x,
          c_x,
          u: u._data ? u._data : u,
        });

        const [deltaMin, j0] = isNegative(delta_x);
        if (j0 === null) {
          executionSteps.push({
            iteration,
            message: "Отрицательных элементов в delta_x не найдено - решение найдено!",
          });
          return { x, J_on, J_ast };
        }

        executionSteps.push({
          iteration,
          message: `Найден отрицательный элемент в delta_x на позиции ${j0}: ${deltaMin}`,
          j0,
        });

        const l = createVectorL(D, B, Ab, J_ast, j0, A);
        const beta = math.dot(math.dot(l, D), math.transpose(l));
        const tetas = getTetas(beta, j0, delta_x, J_ast, l, x);
        const [teta0, teta_j0] = minTetasAndIndex(tetas);

        executionSteps.push({
          iteration,
          message: "Вычислен вектор направления l и значения theta",
          l,
          tetas,
          teta0,
          teta_j0,
        });

        if (teta0 === Infinity) {
          executionSteps.push({
            iteration,
            message: "Задача неограничена - решения не существует",
          });
          return null;
        }

        x = x.map((xi, i) => xi + teta0 * l[i]);
        update(J_on, J_ast, j0, teta_j0, B, A);

        executionSteps.push({
          iteration,
          message: `Обновлены x и множества индексов. Новый x: [${x.map(v => v.toFixed(4)).join(", ")}]`,
          new_J_on: [...J_on],
          new_J_ast: [...J_ast],
        });
      }
    };

    const finalResult = quadraticProgramming(A, D, b, c, J_on, J_ast, x);
    setResult(finalResult);
    setSteps(executionSteps);
  };

  const renderMatrixInput = (matrix, paramName) => (
    <TableContainer component={Paper} sx={{ maxWidth: "fit-content", margin: "10px 0" }}>
      <Table size="small">
        <TableHead>
          <TableRow>
            <TableCell>{paramName}</TableCell>
            {matrix[0].map((_, colIndex) => (
              <TableCell key={`col-${colIndex}`}>Столбец {colIndex + 1}</TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {matrix.map((row, rowIndex) => (
            <TableRow key={`row-${rowIndex}`}>
              <TableCell>Строка {rowIndex + 1}</TableCell>
              {row.map((value, colIndex) => (
                <TableCell key={`cell-${rowIndex}-${colIndex}`}>
                  <TextField
                    type="number"
                    value={value}
                    onChange={(e) => handleParamChange(paramName, rowIndex, colIndex, e.target.value)}
                    size="small"
                    inputProps={{ style: { width: "60px" } }}
                  />
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  );

  const renderVectorInput = (vector, paramName) => (
    <Box sx={{ display: "flex", gap: 1, alignItems: "center", margin: "10px 0" }}>
      <Typography>{paramName}:</Typography>
      {vector.map((value, index) => (
        <TextField
          key={`${paramName}-${index}`}
          type="number"
          value={value}
          onChange={(e) => handleParamChange(paramName, index, undefined, e.target.value)}
          size="small"
          inputProps={{ style: { width: "60px" } }}
        />
      ))}
    </Box>
  );

  const renderIndexSetInput = (set, paramName) => (
    <Box sx={{ display: "flex", gap: 1, alignItems: "center", margin: "10px 0" }}>
      <Typography>{paramName}:</Typography>
      {set.map((value, index) => (
        <TextField
          key={`${paramName}-${index}`}
          type="number"
          value={value}
          onChange={(e) => {
            const newSet = [...set];
            newSet[index] = Number(e.target.value);
            setParams({ ...params, [paramName]: newSet });
          }}
          size="small"
          inputProps={{ style: { width: "60px" } }}
        />
      ))}
      <Button
        variant="outlined"
        size="small"
        onClick={() => {
          const newSet = [...set, 0];
          setParams({ ...params, [paramName]: newSet });
        }}
      >
        Добавить
      </Button>
      {set.length > 0 && (
        <Button
          variant="outlined"
          size="small"
          onClick={() => {
            const newSet = [...set];
            newSet.pop();
            setParams({ ...params, [paramName]: newSet });
          }}
        >
          Удалить
        </Button>
      )}
    </Box>
  );

  return (
    <Container maxWidth="lg" sx={{ my: 4 }}>
      <Paper elevation={3} sx={{ p: 4 }}>
        <Typography variant="h3" component="h1" gutterBottom align="center" sx={{ fontWeight: "bold", mb: 4 }}>
          Решатель задач квадратичного программирования
        </Typography>

        <Typography variant="h5" component="h2" gutterBottom sx={{ mt: 3 }}>
          Постановка задачи
        </Typography>
        <Typography variant="body1" paragraph>
          Задача квадратичного программирования определяется как:
        </Typography>
        <Typography variant="body1" paragraph sx={{ fontFamily: "monospace", textAlign: "center" }}>
          минимизировать (1/2)xᵀDx + cᵀx при условиях Ax = b, x ≥ 0
        </Typography>

        <Divider sx={{ my: 3 }} />

        <Typography variant="h5" component="h2" gutterBottom>
          Входные параметры
        </Typography>

        <Grid container spacing={3}>
          <Grid item xs={12} md={6}>
            <Typography variant="h6">Матрица D (Квадратичный член)</Typography>
            {renderMatrixInput(params.D, "D")}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="h6">Матрица A (Ограничения)</Typography>
            {renderMatrixInput(params.A, "A")}
          </Grid>
          <Grid item xs={12} md={4}>
            <Typography variant="h6">Вектор b (Правая часть ограничений)</Typography>
            {renderVectorInput(params.b, "b")}
          </Grid>
          <Grid item xs={12} md={4}>
            <Typography variant="h6">Вектор c (Линейный член)</Typography>
            {renderVectorInput(params.c, "c")}
          </Grid>
          <Grid item xs={12} md={4}>
            <Typography variant="h6">Начальное x</Typography>
            {renderVectorInput(params.x, "x")}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="h6">Начальное J_on (Активное множество)</Typography>
            {renderIndexSetInput(params.J_on, "J_on")}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="h6">Начальное J_* (Дополнительное множество)</Typography>
            {renderIndexSetInput(params.J_ast, "J_ast")}
          </Grid>
        </Grid>

        <Box sx={{ display: "flex", justifyContent: "center", mt: 4 }}>
          <Button
            variant="contained"
            size="large"
            onClick={handleRun}
            sx={{ px: 6, py: 2, fontSize: "1.1rem" }}
          >
            Решить задачу квадратичного программирования
          </Button>
        </Box>

        {steps.length > 0 && (
          <>
            <Divider sx={{ my: 4 }} />
            <Typography variant="h5" component="h2" gutterBottom>
              Процесс решения
            </Typography>

            <Stepper activeStep={activeStep} orientation="vertical" sx={{ my: 3 }}>
              {steps.map((step, index) => (
                <Step key={`step-${index}`}>
                  <StepLabel onClick={() => setActiveStep(index)} sx={{ cursor: "pointer" }}>
                    {step.message}
                  </StepLabel>
                </Step>
              ))}
            </Stepper>

            {activeStep < steps.length && (
              <Paper elevation={2} sx={{ p: 3, mb: 3 }}>
                <Typography variant="h6" gutterBottom>
                  Итерация {steps[activeStep].iteration}: {steps[activeStep].message}
                </Typography>

                {steps[activeStep].x && (
                  <>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}>Текущий x:</Typography>
                    <Box sx={{ display: "flex", gap: 1 }}>
                      {steps[activeStep].x.map((val, i) => (
                        <Paper key={`x-${i}`} sx={{ p: 1, minWidth: 60, textAlign: "center" }}>
                          x<sub>{i}</sub> = {val.toFixed(4)}
                        </Paper>
                      ))}
                    </Box>
                  </>
                )}

                {steps[activeStep].delta_x && (
                  <>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}>Delta x:</Typography>
                    <Box sx={{ display: "flex", gap: 1 }}>
                      {steps[activeStep].delta_x.map((val, i) => (
                        <Paper 
                          key={`delta-${i}`} 
                          sx={{ 
                            p: 1, 
                            minWidth: 60, 
                            textAlign: "center",
                            backgroundColor: val < 0 ? "#ffebee" : "inherit"
                          }}
                        >
                          Δ<sub>{i}</sub> = {val.toFixed(4)}
                        </Paper>
                      ))}
                    </Box>
                  </>
                )}

                {steps[activeStep].l && (
                  <>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}>Вектор направления l:</Typography>
                    <Box sx={{ display: "flex", gap: 1 }}>
                      {steps[activeStep].l.map((val, i) => (
                        <Paper key={`l-${i}`} sx={{ p: 1, minWidth: 60, textAlign: "center" }}>
                          l<sub>{i}</sub> = {val.toFixed(4)}
                        </Paper>
                      ))}
                    </Box>
                  </>
                )}

                {steps[activeStep].tetas && (
                  <>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}>Шаги (θ):</Typography>
                    <Box sx={{ display: "flex", gap: 1 }}>
                      {steps[activeStep].tetas.map((val, i) => (
                        <Paper 
                          key={`teta-${i}`} 
                          sx={{ 
                            p: 1, 
                            minWidth: 60, 
                            textAlign: "center",
                            backgroundColor: val === steps[activeStep].teta0 ? "#e8f5e9" : "inherit"
                          }}
                        >
                          θ<sub>{i}</sub> = {val === Infinity ? "∞" : val.toFixed(4)}
                        </Paper>
                      ))}
                    </Box>
                  </>
                )}

                {steps[activeStep].new_J_on && (
                  <>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}>Обновленные множества индексов:</Typography>
                    <Box sx={{ display: "flex", gap: 2, mt: 1 }}>
                      <Box>
                        <Typography>J_on: [{steps[activeStep].new_J_on.join(", ")}]</Typography>
                      </Box>
                      <Box>
                        <Typography>J_*: [{steps[activeStep].new_J_ast.join(", ")}]</Typography>
                      </Box>
                    </Box>
                  </>
                )}
              </Paper>
            )}

            <Box sx={{ display: "flex", justifyContent: "space-between", mt: 2 }}>
              <Button
                variant="outlined"
                disabled={activeStep === 0}
                onClick={() => setActiveStep((prev) => prev - 1)}
              >
                Предыдущий шаг
              </Button>
              <Button
                variant="outlined"
                disabled={activeStep >= steps.length - 1}
                onClick={() => setActiveStep((prev) => prev + 1)}
              >
                Следующий шаг
              </Button>
            </Box>
          </>
        )}

        {result && (
          <>
            <Divider sx={{ my: 4 }} />
            <Typography variant="h5" component="h2" gutterBottom>
              Результат
            </Typography>
            <Paper elevation={2} sx={{ p: 3, backgroundColor: "#e8f5e9" }}>
              <Grid container spacing={2}>
                <Grid item xs={12} md={6}>
                  <Typography variant="h6">Оптимальный x:</Typography>
                  <Box sx={{ display: "flex", gap: 1, flexWrap: "wrap", mt: 1 }}>
                    {result.x.map((val, i) => (
                      <Paper key={`result-x-${i}`} sx={{ p: 1, minWidth: 60, textAlign: "center" }}>
                        x<sub>{i}</sub> = {val.toFixed(4)}
                      </Paper>
                    ))}
                  </Box>
                </Grid>
                <Grid item xs={12} md={6}>
                  <Typography variant="h6">Активное множество J_on:</Typography>
                  <Typography sx={{ mt: 1 }}>[{result.J_on.join(", ")}]</Typography>
                </Grid>
                <Grid item xs={12} md={6}>
                  <Typography variant="h6">Дополнтельное множество J_*:</Typography>
                  <Typography sx={{ mt: 1 }}>[{result.J_ast.join(", ")}]</Typography>
                </Grid>
              </Grid>
            </Paper>
          </>
        )}
      </Paper>
    </Container>
  );
};

export default QuadraticProgrammingApp;