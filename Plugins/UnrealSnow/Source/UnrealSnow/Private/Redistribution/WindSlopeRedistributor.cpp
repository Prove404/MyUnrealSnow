#include "Redistribution/WindSlopeRedistributor.h"

using namespace UnrealSnow::Redistribution;

void FWindSlopeRedistributor::EnsureSize(int32 NumCells)
{
	if (WorkA.Num() != NumCells) { WorkA.SetNumZeroed(NumCells); }
	if (WorkB.Num() != NumCells) { WorkB.SetNumZeroed(NumCells); }
	if (GradZ.Num() != NumCells) { GradZ.SetNumZeroed(NumCells); }
	if (GradH.Num() != NumCells) { GradH.SetNumZeroed(NumCells); }
}

void FWindSlopeRedistributor::ComputeGradients(const FSnowGrid& G)
{
	const int32 W = G.W, H = G.H;
	// Central differences for interior, forward/backward for edges.
	for (int32 y = 0; y < H; ++y)
	{
		const int32 y0 = (y > 0) ? (y - 1) : y;
		const int32 y1 = (y < H - 1) ? (y + 1) : y;
		for (int32 x = 0; x < W; ++x)
		{
			const int32 x0 = (x > 0) ? (x - 1) : x;
			const int32 x1 = (x < W - 1) ? (x + 1) : x;
			const int32 i = Idx(x, y, W);
			const float zx0 = G.GroundH_M[Idx(x0, y, W)];
			const float zx1 = G.GroundH_M[Idx(x1, y, W)];
			const float zy0 = G.GroundH_M[Idx(x, y0, W)];
			const float zy1 = G.GroundH_M[Idx(x, y1, W)];
			const float hx0 = G.SnowDepth_M[Idx(x0, y, W)];
			const float hx1 = G.SnowDepth_M[Idx(x1, y, W)];
			const float hy0 = G.SnowDepth_M[Idx(x, y0, W)];
			const float hy1 = G.SnowDepth_M[Idx(x, y1, W)];
			const float inv2dx = 0.5f / FMath::Max(1e-6f, G.CellSizeM);
			GradZ[i].X = (zx1 - zx0) * inv2dx;
			GradZ[i].Y = (zy1 - zy0) * inv2dx;
			GradH[i].X = (hx1 - hx0) * inv2dx;
			GradH[i].Y = (hy1 - hy0) * inv2dx;
		}
	}
}

void FWindSlopeRedistributor::SlopeCreep(FSnowGrid& G, const FRedistributionParams& P, float TanPhi)
{
	const int32 W = G.W, H = G.H;
	const float dt = P.DtHours;
	const float kappa = P.Diffusivity; // m^2/h
	const float maxDelta = P.MaxErodeDepositRate;

	// Total slope of free surface S = grad(Z + H)
	for (int32 i = 0; i < W * H; ++i)
	{
		const FVector2f Slope = FVector2f(GradZ[i].X + GradH[i].X, GradZ[i].Y + GradH[i].Y);
		const float mag = FMath::Sqrt(FMath::Max(0.0f, Slope.X * Slope.X + Slope.Y * Slope.Y));
		WorkA[i] = mag; // store |S| for reuse in wind stage
	}

	// Diffusive creep down the gradient of free-surface, limited when below angle of repose
	// 5-point stencil flux approximation
	for (int32 y = 0; y < H; ++y)
	{
		for (int32 x = 0; x < W; ++x)
		{
			const int32 i = Idx(x, y, W);
			const float slopeMag = WorkA[i];
			float center = G.SnowDepth_M[i];
			float delta = 0.0f;

			// neighbor indices
			const int32 xm = (x > 0) ? Idx(x - 1, y, W) : i;
			const int32 xp = (x < W - 1) ? Idx(x + 1, y, W) : i;
			const int32 ym = (y > 0) ? Idx(x, y - 1, W) : i;
			const int32 yp = (y < H - 1) ? Idx(x, y + 1, W) : i;

			// Laplacian of (Z+H) approximated via H only for redistribution amount
			const float lap = (G.SnowDepth_M[xm] + G.SnowDepth_M[xp] + G.SnowDepth_M[ym] + G.SnowDepth_M[yp] - 4.0f * center);
			const float diffused = kappa * dt * lap / FMath::Max(1e-6f, G.CellSizeM * G.CellSizeM);
			// If above angle of repose, enhance downslope transport proportional to excess
			const float tanRatio = (TanPhi > 1e-6f) ? FMath::Min(1.0f, slopeMag / TanPhi) : 1.0f;
			delta = diffused * tanRatio;
			// Clamp per-iteration change
			WorkB[i] = FMath::Clamp(delta, -maxDelta, maxDelta);
		}
	}

	// Apply with conservation: accumulate neighbor exchanges explicitly
	// Here we perform a simple Jacobi step where WorkB stores net change at cell
	for (int32 i = 0; i < W * H; ++i)
	{
		G.SnowDepth_M[i] = FMath::Max(0.0f, G.SnowDepth_M[i] + WorkB[i]);
	}
}

void FWindSlopeRedistributor::WindAdvection(FSnowGrid& G, const FRedistributionParams& P, float TanPhi)
{
	const int32 W = G.W, H = G.H;
	const float dt = P.DtHours;
	const float speed = P.WindSpeedMS;
	FVector2f wdir = FVector2f((float)P.WindDir.X, (float)P.WindDir.Y);
	const float wmag = FMath::Sqrt(FMath::Max(1e-8f, wdir.X * wdir.X + wdir.Y * wdir.Y));
	wdir.X /= wmag; wdir.Y /= wmag;

	// Capacity limited transport
	for (int32 y = 0; y < H; ++y)
	{
		for (int32 x = 0; x < W; ++x)
		{
			const int32 i = Idx(x, y, W);
			const FVector2f nSlopeZ = [] (const FVector2f& g)
			{
				const float m = FMath::Sqrt(FMath::Max(1e-8f, g.X * g.X + g.Y * g.Y));
				return FVector2f(g.X / m, g.Y / m);
			}(GradZ[i]);

			const float dotW = nSlopeZ.X * wdir.X + nSlopeZ.Y * wdir.Y;
			const float Exposure = FMath::Clamp(-dotW, 0.0f, 1.0f); // lee-side shelter
			const float SpeedUp = FMath::Clamp(dotW, 0.0f, 1.0f);   // windward acceleration

			const float SMag = FMath::Sqrt(FMath::Max(0.0f, (GradZ[i].X + GradH[i].X) * (GradZ[i].X + GradH[i].X) + (GradZ[i].Y + GradH[i].Y) * (GradZ[i].Y + GradH[i].Y)));
			const float angleFactor = 1.0f - FMath::Min(1.0f, (TanPhi > 1e-6f) ? (SMag / TanPhi) : 0.0f);

			const float speedExcess = FMath::Max(0.0f, speed - P.UstarThresh);
			const float Q = P.CapacityK * FMath::Pow(speedExcess, P.CapacityP) * angleFactor; // [m/h]
			const float move = FMath::Min(P.MaxErodeDepositRate, Q * dt);
			const float available = FMath::Max(0.0f, G.SnowDepth_M[i]);
			const float transported = FMath::Min(move, available);

			// Move downwind by 1 cell via bilinear weights to 2 neighbor cells along primary axis
			const float fx = wdir.X; const float fy = wdir.Y;
			const int32 tx = FMath::Clamp((fx >= 0.0f) ? (x + 1) : (x - 1), 0, W - 1);
			const int32 ty = FMath::Clamp((fy >= 0.0f) ? (y + 1) : (y - 1), 0, H - 1);
			const float ax = FMath::Abs(fx); const float ay = FMath::Abs(fy);

			WorkB[i] = -transported * (0.5f + 0.5f * SpeedUp); // erosion weight on windward
			// distribute to two targets proportional to direction components
			const float toX = transported * ax * (0.5f + 0.5f * Exposure);
			const float toY = transported * ay * (0.5f + 0.5f * Exposure);
			WorkA[Idx(tx, y, W)] += toX;
			WorkA[Idx(x, ty, W)] += toY;
		}
	}

	// Apply erosion/deposition and clamp
	for (int32 i = 0; i < W * H; ++i)
	{
		G.SnowDepth_M[i] = FMath::Max(0.0f, G.SnowDepth_M[i] + WorkB[i] + WorkA[i]);
		WorkA[i] = 0.0f; // reset for reuse next substep
		WorkB[i] = 0.0f;
	}
}

void FWindSlopeRedistributor::Step(FSnowGrid& Grid, const FRedistributionParams& P)
{
	if (Grid.W <= 1 || Grid.H <= 1) return;
	const int32 Num = Grid.W * Grid.H;
	if (Grid.GroundH_M.Num() != Num || Grid.SnowDepth_M.Num() != Num) return;

	EnsureSize(Num);

	const float TanPhi = FMath::Tan(FMath::DegreesToRadians(P.AngleOfReposeDeg));
	for (int32 iter = 0; iter < FMath::Max(1, P.Iterations); ++iter)
	{
		ComputeGradients(Grid);
		SlopeCreep(Grid, P, TanPhi);
		WindAdvection(Grid, P, TanPhi);
	}
}


