# [Model Components](@id model_components)

The model is structured using model components as composite types. Here is a general overview of the model components:

```mermaid
graph TD
  A[ParameterSet] --> B[SurfaceFractions] --> F[LocationSpecificSurfaceFractions]
  A --> C[ThermalProperties]
  A --> D[UrbanGeometryParameters]
  A --> E[VegetationParameters] --> G[HeightDependentVegetationParameters]
```
