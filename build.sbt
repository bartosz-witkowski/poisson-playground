val scala3Version = "3.2.0"

ThisBuild / run / javaOptions += "--add-modules=jdk.incubator.vector"

scalacOptions ++= List(
  "-Xjline:off")

lazy val root = project
  .in(file("."))
  .settings(
    name := "poisson",
    version := "0.1.0-SNAPSHOT",
    scalaVersion := scala3Version,
    libraryDependencies ++= List(
        "org.scala-lang.modules" %% "scala-parallel-collections" % "1.0.4",
	"org.scalameta" %% "munit" % "0.7.29" % Test)
  )
